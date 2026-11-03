"""Source-only closed-box receptor coverage and separately frozen 8 A repair.

Legacy input archives are read-only. No pose, score, RMSD or recovery file is
opened by audit/preparation. Corrected docking is a separate explicit command.
"""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import re
from pathlib import Path
import shutil
import tempfile

import gemmi
import numpy as np

from .aav import PROTEIN_RESIDUES, _pdb_atoms, prepare_receptor
from .provenance import sha256, write_json, write_text_atomic
from .site_geometry import read_receptor
from .site_prediction_verification import restore_site_inputs

LEGACY_CONFIG = Path('benchmarks/blind_pose_experiment.json')
AMENDMENT = Path('benchmarks/blind_pose_neighborhood_amendment.json')
AUDIT = Path('results/blind_pose/neighborhood_coverage_audit.json')
CORRECTED_RESULTS = Path('results/blind_pose_neighborhood8')
CORRECTED_PREP = Path('data/processed/blind_pose_neighborhood8')
EXPECTED_BOXES = 52
STANDARD_CUTOFF = 8.0


def _load(path):
    return json.loads(Path(path).read_text())


def _raw_sha(raw):
    return hashlib.sha256(raw).hexdigest()


def closed_box_distances(coordinates, center, size):
    """Vectorized Euclidean atom distance to a closed axis-aligned search box."""
    xyz = np.asarray(coordinates, dtype=float)
    center, size = np.asarray(center, dtype=float), np.asarray(size, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1:] != (3,) or center.shape != (3,) or size.shape != (3,):
        raise ValueError('Coordinates require Nx3 and box requires two three-vectors')
    if not np.isfinite(xyz).all() or not np.isfinite(center).all() or not np.isfinite(size).all() or np.any(size <= 0):
        raise ValueError('Finite coordinates and positive box dimensions required')
    return np.linalg.norm(np.maximum(np.abs(xyz - center) - size / 2, 0.), axis=1)


def chain_selection(keys, coordinates, center, size, cutoff):
    """Retain whole chains if any of their source atoms is within cutoff."""
    if not np.isfinite(cutoff) or cutoff <= 0 or len(keys) != len(coordinates) or not keys:
        raise ValueError('Nonempty atom inventory and positive finite cutoff required')
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate source heavy-atom identity')
    distances = closed_box_distances(coordinates, center, size)
    chains = np.asarray([k[0] for k in keys])
    unique, inverse = np.unique(chains, return_inverse=True)
    minima = np.full(len(unique), np.inf)
    np.minimum.at(minima, inverse, distances)
    selected = unique[minima <= cutoff].tolist()
    if not selected:
        raise ValueError('No source chains intersect the expanded box')
    mask = np.isin(chains, selected)
    omitted = np.flatnonzero(~mask)
    nearest_index = int(omitted[np.argmin(distances[omitted])]) if len(omitted) else None
    witness = None if nearest_index is None else dict(atom_key=list(keys[nearest_index]),
            coordinates_A=np.asarray(coordinates[nearest_index]).tolist(), distance_to_closed_box_A=float(distances[nearest_index]))
    return dict(cutoff_A=float(cutoff), selected_source_chains=selected,
                selected_heavy_atoms=int(mask.sum()), source_heavy_atoms=len(keys),
                source_chain_count=len(unique), nearest_omitted_atom=witness,
                closest_omitted_atom_to_box_A=None if witness is None else witness['distance_to_closed_box_A'],
                minimum_distance_per_source_chain_A={str(k): float(v) for k, v in zip(unique, minima)})


def _archive_inputs(directory):
    """Verify and decode only declared preparation-input bytes, never outcomes."""
    directory = Path(directory)
    manifest_path = directory / 'manifest.json'
    manifest = _load(manifest_path)
    content, receipts = {}, {}
    for item in manifest['files']:
        if any(Path(item[k]).name != item[k] for k in ['filename', 'archive']):
            raise ValueError('Input archive path traversal')
        source = directory / item['archive']
        if source.is_symlink() or sha256(source) != item['archive_sha256']:
            raise ValueError('Prepared input compressed checksum mismatch')
        raw = gzip.decompress(source.read_bytes())
        if _raw_sha(raw) != item['sha256'] or len(raw) != item['original_bytes']:
            raise ValueError('Prepared input byte checksum mismatch')
        content[item['filename']] = raw
        receipts[item['filename']] = dict(sha256=item['sha256'], archive_sha256=item['archive_sha256'], bytes=len(raw))
    if len(content) != len(manifest['files']):
        raise ValueError('Duplicate input archive filename')
    return content, dict(manifest_sha256=sha256(manifest_path), files=receipts)


def _heavy_proof(keys, coordinates, pdb_text, pdbqt_text, selected):
    source = dict(zip(keys, coordinates))
    expected = {k for k in keys if k[0] in selected}
    results = {}
    for label, text, qt in [('PDB', pdb_text, False), ('PDBQT', pdbqt_text, True)]:
        parsed = {k: a for k, a in _pdb_atoms(text, pdbqt=qt).items() if not a['hydrogen']}
        if set(parsed) != expected:
            raise ValueError(f'{label} source heavy-atom identities differ from complete selected chains')
        displacement = np.asarray([np.asarray(parsed[k]['coordinates']) - source[k] for k in sorted(expected)])
        maximum = float(np.linalg.norm(displacement, axis=1).max())
        if maximum > .002:
            raise ValueError(f'{label} source heavy coordinates moved beyond existing preparation tolerance')
        serialised = [[*k, *parsed[k]['coordinates']] for k in sorted(parsed)]
        results[label] = dict(heavy_atoms=len(parsed), source_heavy_identity_exact=True,
                             source_heavy_coordinates_exact=bool(maximum == 0), max_coordinate_displacement_A=maximum,
                             identity_coordinate_sha256=_raw_sha(json.dumps(serialised, separators=(',', ':')).encode()))
    return results


def _source_metadata(root, condition, restored):
    data = _load(restored / 'receptor.cif.manifest.json')
    if data['source_sha256'] != condition['source_sha256'] or data['receptor_sha256'] != sha256(restored / 'receptor.cif'):
        raise ValueError('Original site-source provenance differs from frozen docking source')
    return dict(original_source_path=condition['source'], original_source_sha256=condition['source_sha256'],
                source_receptor_sha256=data['receptor_sha256'],
                original_site_input_manifest=f'results/site_discovery/{condition["condition_id"]}/inputs/manifest.json',
                original_site_input_manifest_sha256=sha256(root / f'results/site_discovery/{condition["condition_id"]}/inputs/manifest.json'),
                stripping_policy=data['policy'], assembly=data['assembly'])


def audit_legacy(root):
    root = Path(root)
    config = _load(root / LEGACY_CONFIG)
    config_sha = sha256(root / LEGACY_CONFIG)
    freeze = _load(root / 'results/blind_pose/freeze_receipt.json')
    if freeze['configuration_sha256'] != config_sha or config['receptor_policy']['interaction_cutoff_A'] != 6.:
        raise ValueError('Legacy frozen six-angstrom configuration differs')
    paths = sorted((root / 'results/blind_pose').glob('*/*/box*/inputs/manifest.json'))
    if len(paths) != EXPECTED_BOXES:
        raise ValueError('Legacy prepared-box denominator differs from 52')
    by_condition = {c['condition_id']: c for c in config['cases']}
    rows, condition_preparations, expected = [], [], set()
    for condition in config['cases']:
        for method in config['methods']:
            path = root / 'results/blind_pose' / condition['condition_id'] / method / 'condition_preparation.json'
            preparation = _load(path)
            jobs = [(j['id'], j['status']) for j in preparation['jobs']]
            expected.update((condition['condition_id'], method, name) for name, state in jobs if state == 'prepared')
            condition_preparations.append(dict(condition_id=condition['condition_id'], method=method,
                status=preparation['status'], jobs=[dict(box_id=n, preparation_status=s) for n, s in jobs],
                preparation_manifest_sha256=sha256(path)))
    observed = {(p.parents[3].name, p.parents[2].name, p.parents[1].name) for p in paths}
    if expected != observed:
        raise ValueError('Prepared condition manifests and archived 52 boxes disagree')
    with tempfile.TemporaryDirectory(prefix='neighborhood-source-audit-') as tmp:
        cache = {}
        for path in paths:
            condition_id, method, box_id = path.parents[3].name, path.parents[2].name, path.parents[1].name
            if condition_id not in cache:
                dest = Path(tmp) / condition_id
                restore_site_inputs(root, condition_id, 'original', dest)
                keys, xyz, elements = read_receptor(dest / 'receptor.cif')
                if any(k[3] not in PROTEIN_RESIDUES for k in keys):
                    raise ValueError('Nonprotein atoms in source discovery receptor')
                cache[condition_id] = (keys, xyz, _source_metadata(root, by_condition[condition_id], dest))
            keys, xyz, source = cache[condition_id]
            raw, receipts = _archive_inputs(path.parent)
            prep = json.loads(raw['preparation.json'])
            if prep['configuration_sha256'] != config_sha or prep['source_receptor_sha256'] != source['source_receptor_sha256']:
                raise ValueError('Prepared input does not bind frozen source/configuration')
            six = chain_selection(keys, xyz, prep['box_center'], prep['box_size'], 6.)
            eight = chain_selection(keys, xyz, prep['box_center'], prep['box_size'], 8.)
            if six['selected_source_chains'] != prep['neighborhood']['selected_source_chains']:
                raise ValueError('Legacy source chain selection does not reproduce frozen six-angstrom rule')
            old_nearest = prep['neighborhood']['closest_omitted_atom_to_box_A']
            new_nearest = six['closest_omitted_atom_to_box_A']
            if (old_nearest is None) != (new_nearest is None) or (old_nearest is not None and abs(old_nearest - new_nearest) > 1e-10):
                raise ValueError('Legacy omitted-distance witness differs')
            proof = _heavy_proof(keys, xyz, raw['neighborhood.pdb'].decode(), raw['neighborhood.pdbqt'].decode(), six['selected_source_chains'])
            if _raw_sha(raw['ligand.pdbqt']) != prep['ligand_sha256']:
                raise ValueError('Independent ligand input changed')
            atom_types = sorted({line[77:].strip() for line in raw['ligand.pdbqt'].decode().splitlines() if line.startswith(('ATOM  ', 'HETATM'))})
            if any(re.match(r'^(?:G|CG)\d', t) for t in atom_types):
                raise ValueError('Macrocycle glue atom type requires separate cutoff audit')
            identity = f'{condition_id}/{method}/{box_id}'
            added = sorted(set(eight['selected_source_chains']) - set(six['selected_source_chains']))
            rows.append(dict(box_key=identity, condition_id=condition_id, complex_id=prep['complex_id'], method=method,
                box_id=box_id, box_center_A=prep['box_center'], box_size_A=prep['box_size'],
                legacy_configuration_sha256=config_sha, source=source, archived_inputs=receipts,
                source_heavy_atom_proof=proof, original_ligand_sha256=prep['ligand_sha256'], ligand_Vina_atom_types=atom_types,
                selected_at_6A=six, selected_at_8A=eight, added_chains_for_8A=added,
                added_heavy_atoms_for_8A=eight['selected_heavy_atoms']-six['selected_heavy_atoms'],
                legacy_6A_input_covers_standard_Vina_8A=not added,
                whole_source_receptor_retained=six['source_heavy_atoms']==six['selected_heavy_atoms'],
                native_coordinate_access_for_this_coverage_audit=False,
                center_information_regime='frozen_oracle_center' if method == 'oracle' else 'frozen_predictor_only_center'))
    affected = [r['box_key'] for r in rows if not r['legacy_6A_input_covers_standard_Vina_8A']]
    return dict(schema_version=1, source_only=True, docking_or_accuracy_outcomes_consulted=False,
        algorithm='vectorized Euclidean distance to closed box; whole-chain minima <= cutoff retained',
        algorithm_code_sha256=sha256(Path(__file__)), legacy_configuration=str(LEGACY_CONFIG),
        legacy_configuration_sha256=config_sha, prepared_box_count=len(rows), condition_method_count=len(condition_preparations),
        AAV_prepared_boxes=sum(r['complex_id']=='aav9' for r in rows),
        external_prepared_boxes=sum(r['complex_id']!='aav9' for r in rows),
        legacy_boxes_covering_8A=sum(r['legacy_6A_input_covers_standard_Vina_8A'] for r in rows),
        affected_box_keys=affected, scope='Development geometry repair. Coverage is relative to the archived dry protein-only source receptor, not certification of biological assembly completeness.',
        preparation_conditions=condition_preparations, boxes=rows)


def _verify_amendment(root):
    amendment = _load(root / AMENDMENT)
    if amendment['legacy_configuration_sha256'] != sha256(root / LEGACY_CONFIG):
        raise ValueError('Legacy configuration changed')
    if amendment['legacy_coverage_audit_sha256'] != sha256(root / AUDIT):
        raise ValueError('Source-only coverage audit changed after amendment freeze')
    for item in amendment['authoritative_cutoff_sources']:
        if sha256(root / item['archive']) != item['sha256']:
            raise ValueError('Authoritative cutoff source changed')
    software = _load(root / CORRECTED_RESULTS / 'software/adapter_manifest.json')
    package = Path(__file__).parent
    for name, digest in software['package_files'].items():
        if sha256(package / name) != digest:
            raise ValueError('Run preparation through frozen software snapshot; current helper changed: ' + name)
    return amendment


def freeze_amendment(root):
    """Freeze source reason and exact correction denominator before new runs."""
    root = Path(root)
    if (root / AMENDMENT).exists():
        raise FileExistsError('Preserve existing amendment freeze')
    audit = _load(root / AUDIT)
    if audit['algorithm_code_sha256'] != sha256(Path(__file__)) or audit['prepared_box_count'] != EXPECTED_BOXES:
        raise ValueError('Coverage audit must be current before freezing')
    software = root / CORRECTED_RESULTS / 'software'
    sources = _load(software / 'cutoff_sources.json')['sources']
    package = software / 'aav_benchmark';package.mkdir(parents=True, exist_ok=False)
    package_files = {}
    for path in sorted(Path(__file__).parent.glob('*.py')):
        shutil.copyfile(path, package / path.name);package_files[path.name] = sha256(path)
    scriptdir = software / 'scripts';scriptdir.mkdir()
    driver = root / 'scripts/audit_neighborhood_coverage.py'
    shutil.copyfile(driver, scriptdir / driver.name)
    write_json(software / 'adapter_manifest.json', dict(package_files=package_files,
        driver_sha256=sha256(driver), policy='Pinned complete local package and separate driver; no monkeypatch of legacy helpers'))
    amendment = dict(schema_version=1, freeze_date='2026-11-03',
        scope='Development receptor-neighborhood geometry correction and paired sensitivity; not an untouched confirmatory study',
        frozen_before_corrected_docking=True, new_AAV_accuracy_outcomes_consulted=False,
        legacy_configuration=str(LEGACY_CONFIG), legacy_configuration_sha256=sha256(root / LEGACY_CONFIG),
        legacy_coverage_audit=str(AUDIT), legacy_coverage_audit_sha256=sha256(root / AUDIT),
        authoritative_cutoff_sources=sources, legacy_cutoff_A=6., corrected_standard_nonbonded_cutoff_A=8.,
        rule='Retain each complete archived source protein chain with any heavy atom at distance <=8 A from the closed search box',
        correction_selection='Every legacy prepared box whose 6 A selected-chain set differs from its 8 A set; source geometry only',
        affected_box_keys=audit['affected_box_keys'], corrected_box_count=len(audit['affected_box_keys']),
        planned_seed_count=3, planned_run_count=3*len(audit['affected_box_keys']),
        preserved_settings=['independently embedded ligand bytes','box center','box dimensions','candidate ranks','protein-only dry rigid model','seeds104729/130363/155921','cpu2','exhaustiveness16','num_modes20','min_rmsd0.5','energy_range10'],
        changed_setting='whole-chain coverage margin only, 6 to8 A; no source heavy-atom repair',
        legacy_preservation='Existing configuration, inputs, raw runs and outcomes remain unchanged; corrections use separate processed/result trees',
        primary_atomic_external_effect='All40 prepared non-AAV boxes retain whole source receptors and their chain sets do not change; this repair cannot alter those inputs',
        AAV_atomic_scope='AAV cases have no qualified atomic reference in this experiment; correction does not earn one',
        untouched_future_workflow='Use standard8 A coverage before first freeze of any new untouched workflow',
        execution_gate='Preparation/preflight allowed now. Start corrected docking only after coordinator grants a CPU slot; no launch is performed by freeze or prepare commands.')
    write_json(root / AMENDMENT, amendment)
    config = copy.deepcopy(_load(root / LEGACY_CONFIG));config['receptor_policy']['interaction_cutoff_A']=8.
    config['receptor_policy']['selection']=amendment['rule'];config['purpose']='Development sensitivity for three source-diagnosed neighborhood defects; preserve all legacy evidence'
    config['cohort_role']='development geometry repair, selected before corrected docking from source coverage'
    config['neighborhood_amendment_sha256']=sha256(root / AMENDMENT)
    write_json(root / CORRECTED_RESULTS / 'configuration.json', config)
    write_json(root / CORRECTED_RESULTS / 'freeze_receipt.json', dict(configuration_sha256=sha256(root / CORRECTED_RESULTS / 'configuration.json'),
        amendment_sha256=sha256(root / AMENDMENT), frozen_before_new_docking=True, new_pose_outcomes_consulted=False,
        corrected_box_keys=amendment['affected_box_keys'], expected_runs=amendment['planned_run_count']))
    return amendment


def _eight_protein_text(source_path, keys, xyz, center, size):
    selected = chain_selection(keys, xyz, center, size, 8.)
    structure = gemmi.read_structure(str(source_path))
    out = gemmi.Structure();out.add_model(gemmi.Model('1'))
    for chain in structure[0]:
        if chain.name in selected['selected_source_chains']:
            if len(chain.name) != 1:
                raise ValueError('Explicit chain map required for PDB serialization')
            out[0].add_chain(chain.clone())
    if selected['selected_heavy_atoms'] > 99999:
        raise ValueError('PDB atom-serial capacity exceeded')
    return out.make_pdb_string(), selected


def _archive_box(prepared, destination):
    destination.mkdir(parents=True, exist_ok=False)
    names=['ligand.pdbqt','neighborhood.pdb','neighborhood.pdbqt','neighborhood_atom_changes.json','preparation.json','preparation_attempt.json','reference.json']
    entries=[]
    for name in names:
        raw=(prepared/name).read_bytes();archive=name+'.gz';target=destination/archive;target.write_bytes(gzip.compress(raw,mtime=0))
        entries.append(dict(filename=name,archive=archive,sha256=_raw_sha(raw),archive_sha256=sha256(target),original_bytes=len(raw)))
    write_json(destination/'manifest.json',dict(schema_version=1,files=entries))
    return sha256(destination/'manifest.json')


def prepare_corrections(root):
    root=Path(root);amendment=_verify_amendment(root);audit=_load(root/AUDIT)
    rows={x['box_key']:x for x in audit['boxes']};config_path=root/CORRECTED_RESULTS/'configuration.json'
    if (root/CORRECTED_RESULTS/'prepared_corrections.json').exists():
        raise FileExistsError('Corrected inputs already frozen; verify rather than overwrite')
    corrections=[]
    with tempfile.TemporaryDirectory(prefix='neighborhood8-prepare-') as tmp:
        cache={}
        for key in amendment['affected_box_keys']:
            row=rows[key];cond,method,box=key.split('/')
            if row['complex_id']!='aav9':
                raise ValueError('This frozen correction scope contains AAV boxes only')
            if cond not in cache:
                source=Path(tmp)/cond;restore_site_inputs(root,cond,'original',source)
                keys,xyz,elements=read_receptor(source/'receptor.cif');cache[cond]=(source,keys,xyz)
            source,keys,xyz=cache[cond]
            legacy_dir=root/'results/blind_pose'/cond/method/box/'inputs'
            raw,receipt=_archive_inputs(legacy_dir);legacy=json.loads(raw['preparation.json'])
            output=root/CORRECTED_PREP/cond/method/box;output.mkdir(parents=True,exist_ok=False)
            protein,selection=_eight_protein_text(source/'receptor.cif',keys,xyz,row['box_center_A'],row['box_size_A'])
            log=prepare_receptor(protein,output/'neighborhood')
            for name in ['ligand.pdbqt','reference.json']:(output/name).write_bytes(raw[name])
            prep=copy.deepcopy(legacy);prep.pop('receptor_cache',None)
            prep.update(configuration_sha256=sha256(config_path),preparation_code_sha256=sha256(Path(__file__)),
                        receptor_preparation=log,legacy_preparation_sha256=_raw_sha(raw['preparation.json']),
                        neighborhood_amendment_sha256=sha256(root/AMENDMENT),
                        legacy_box_key=key,source_receptor_sha256=sha256(source/'receptor.cif'),
                        neighborhood=dict(rule=amendment['rule'],selected_source_chains=selection['selected_source_chains'],
                            source_atom_count=len(keys),source_chain_count=selection['source_chain_count'],
                            selected_heavy_atoms=selection['selected_heavy_atoms'],interaction_cutoff_A=8.,
                            closest_omitted_atom_to_box_A=selection['closest_omitted_atom_to_box_A'],
                            source_heavy_atoms_unchanged=True,native_coordinate_access=False),
                        cohort_role='development_geometry_repair')
            prep['receptors'][0]['assembly_condition']='complete_box_neighborhood8A'
            proof=_heavy_proof(keys,xyz,(output/'neighborhood.pdb').read_text(),(output/'neighborhood.pdbqt').read_text(),selection['selected_source_chains'])
            if sha256(output/'ligand.pdbqt')!=row['original_ligand_sha256']:
                raise ValueError('Correction changed frozen independent ligand bytes')
            if selection['closest_omitted_atom_to_box_A'] is not None and selection['closest_omitted_atom_to_box_A']<=8:
                raise ValueError('Corrected neighborhood remains incomplete at8 A')
            write_json(output/'preparation.json',prep)
            write_json(output/'preparation_attempt.json',dict(status='ok',source_only=True,legacy_box_key=key,
                configuration_sha256=sha256(config_path),amendment_sha256=sha256(root/AMENDMENT)))
            archived=root/CORRECTED_RESULTS/cond/method/box/'inputs';manifest_sha=_archive_box(output,archived)
            corrections.append(dict(box_key=key,prepared_directory=str(output.relative_to(root)),
                preparation_sha256=sha256(output/'preparation.json'),input_archive=str(archived.relative_to(root)),
                input_manifest_sha256=manifest_sha,legacy_input_manifest_sha256=receipt['manifest_sha256'],
                original_ligand_sha256=row['original_ligand_sha256'],selection_8A=selection,
                added_chains=row['added_chains_for_8A'],source_heavy_atom_proof=proof,
                source_receptor_sha256=prep['source_receptor_sha256']))
    manifest=dict(schema_version=1,status='prepared_waiting_for_coordinator_CPU_slot',source_only=True,
        amendment_sha256=sha256(root/AMENDMENT),configuration_sha256=sha256(config_path),
        frozen_before_corrected_docking=True,corrected_box_count=len(corrections),expected_runs=3*len(corrections),
        corrections=corrections)
    write_json(root/CORRECTED_RESULTS/'prepared_corrections.json',manifest)
    return manifest


def verify_corrections(root):
    root=Path(root);amendment=_verify_amendment(root);manifest=_load(root/CORRECTED_RESULTS/'prepared_corrections.json')
    if manifest['amendment_sha256']!=sha256(root/AMENDMENT) or manifest['configuration_sha256']!=sha256(root/CORRECTED_RESULTS/'configuration.json'):
        raise ValueError('Corrected preparation configuration freeze mismatch')
    if [x['box_key'] for x in manifest['corrections']]!=amendment['affected_box_keys']:
        raise ValueError('Corrected box denominator changed')
    with tempfile.TemporaryDirectory(prefix='neighborhood8-replay-') as tmp:
        cache={}
        for row in manifest['corrections']:
            cond,method,box=row['box_key'].split('/')
            if cond not in cache:
                dest=Path(tmp)/cond;restore_site_inputs(root,cond,'original',dest)
                keys,xyz,elements=read_receptor(dest/'receptor.cif');cache[cond]=(keys,xyz)
            keys,xyz=cache[cond];raw,receipt=_archive_inputs(root/row['input_archive'])
            if receipt['manifest_sha256']!=row['input_manifest_sha256']:
                raise ValueError('Corrected input archive manifest changed')
            prep=json.loads(raw['preparation.json'])
            if _raw_sha(raw['preparation.json'])!=row['preparation_sha256'] or _raw_sha(raw['ligand.pdbqt'])!=row['original_ligand_sha256']:
                raise ValueError('Corrected preparation or unchanged ligand bytes differ')
            selected=chain_selection(keys,xyz,prep['box_center'],prep['box_size'],8.)
            if selected!=row['selection_8A']:
                raise ValueError('Corrected source geometry replay differs')
            proof=_heavy_proof(keys,xyz,raw['neighborhood.pdb'].decode(),raw['neighborhood.pdbqt'].decode(),selected['selected_source_chains'])
            if proof!=row['source_heavy_atom_proof']:
                raise ValueError('Corrected heavy-atom identity/coordinates differ')
    return dict(status='verified',source_only=True,docking_or_accuracy_outcomes_consulted=False,
                corrected_boxes=len(manifest['corrections']),expected_runs=manifest['expected_runs'])


def run_corrections(root,executable,*,resume=False):
    """Explicit execution entry point; coordinator CPU-slot grant is required."""
    from .engine import dock
    root=Path(root);verify_corrections(root);manifest=_load(root/CORRECTED_RESULTS/'prepared_corrections.json')
    freeze=_load(root/CORRECTED_RESULTS/'freeze_receipt.json');config=root/CORRECTED_RESULTS/'configuration.json'
    lock=_load(root/'benchmarks/vina.lock.json')
    if sha256(executable) not in {x['sha256'] for x in lock['assets'].values()}:
        raise ValueError('Executable is not a pinned official Vina build')
    if freeze['configuration_sha256']!=sha256(config):raise ValueError('Corrected run configuration changed')
    runs=[]
    for correction in manifest['corrections']:
        cond,method,box=correction['box_key'].split('/')
        result=dock(root/correction['prepared_directory'],root/CORRECTED_RESULTS/cond/method/box,executable,config,resume=resume)
        runs.append(dict(box_key=correction['box_key'],runs=result))
        write_json(root/CORRECTED_RESULTS/'execution_manifest.json',dict(schema_version=1,amendment_sha256=sha256(root/AMENDMENT),
            configuration_sha256=sha256(config),jobs=runs,status='ok' if len(runs)==3 and all(x['status']=='ok' for j in runs for x in j['runs']) else 'running_or_failed'))
    return runs
