"""Synthetic TEMP byte packaging only; no real structural inventory or readers."""
from copy import deepcopy
import gzip
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
from unittest.mock import patch

import pytest

SOURCE = Path(__file__).resolve().parents[1] / 'results/structural_pose_output_foundation/package_outputs.py'
SPEC = importlib.util.spec_from_file_location('synthetic_structural_packager', SOURCE)
pkg = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pkg)


@pytest.fixture
def tmp_path(tmp_path_factory):
    # The production API rejects every supplied alias, including system TEMP
    # aliases. Ordinary test roots therefore use their canonical paths.
    return tmp_path_factory.mktemp('structural-output-packager').resolve()


def put(root, name, raw):
    path = root / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(raw)
    return {'path': name, 'bytes': len(raw), 'sha256': pkg.digest(raw)}


def write_json(root, name, value):
    return put(root, name, pkg.encoded(value))


def historical_tar(rows, *, corrupt=False):
    output = io.BytesIO()
    with gzip.GzipFile(fileobj=output, filename='', mode='wb', mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode='w', format=tarfile.PAX_FORMAT) as archive:
            for index, (name, raw) in enumerate(sorted(rows.items())):
                if corrupt and index == 0:
                    raw = b'changed opaque historical bytes'
                info = tarfile.TarInfo(name); info.size = len(raw); info.mtime = 0
                archive.addfile(info, io.BytesIO(raw))
    return output.getvalue()


def fixture(root):
    for name in pkg.SOURCE_ROOTS.values():
        (root / name).mkdir(parents=True)
    put(root, pkg.SOURCE_ROOTS['main'] + '/history/partial.bin', b'opaque\x00partial\xffhistory\n')
    put(root, pkg.SOURCE_ROOTS['main'] + '/empty', b'')
    put(root, pkg.SOURCE_ROOTS['assembly'] + '/failed.stdout', b'opaque failure evidence; no invented telemetry\n')
    prefix = pkg.SOURCE_ROOTS['haddock']; members = {'native/old-a.bin': b'preserved\x00first', 'native/old-b.bin': b'opaque unparsed energy NaN\n'}
    rows = [put(root, prefix + '/' + name, raw) for name, raw in sorted(members.items())]
    manifest = write_json(root, prefix + '/orchestration_recovery/historical_manifest.json', {'schema_version': 1, 'files': rows})
    archive = put(root, prefix + '/orchestration_recovery/historical_raw.tar.gz', historical_tar(members))
    contract = put(root, 'benchmarks/source-protocol.json', b'opaque published protocol\n')
    return {'schema_version': 1, 'profile': 'synthetic_test_fixture', 'source_roots': deepcopy(pkg.SOURCE_ROOTS),
            'destination': 'packaged', 'limits': {'max_raw_bytes': 100_000, 'max_files': 3, 'max_gzip_bytes': 1_000_000},
            'exclusions': [], 'historical_first_job': {'archive': archive, 'manifest': manifest, 'member_count': len(rows), 'target_prefix': prefix},
            'source_contracts': [contract], 'source_protocol_commit': 'a' * 40, 'must_be_absent': list(pkg.MUST_BE_ABSENT)}


def approve(root, plan, *, readiness_changes=None, approval_changes=None):
    plan_pin = write_json(root, 'plan.json', plan)
    ready = {'schema_version': 1, 'status': 'all_structural_sources_terminal_metadata_only', 'accuracy_computed': False,
             'numeric_pose_or_energy_payloads_interpreted': False, 'writers_stopped': True, 'pending': [],
             'counts': dict(pkg.READINESS_COUNTS, retained=2, recorded=2),
             'source_snapshot_sha256': pkg.snapshot_sha256(pkg.inventory(root, plan))}
    ready.update(readiness_changes or {})
    validation_pin = write_json(root, 'readiness.json', ready)
    approval = {'schema_version': 1, 'status': 'approved_terminal_structural_raw_packaging', 'packaging_authorized': True,
                'packaging_plan': plan_pin, 'packager_sha256': pkg.digest(SOURCE.read_bytes()), 'terminal_validation': validation_pin}
    approval.update(approval_changes or {})
    write_json(root, 'approval.json', approval)
    return ready


def run(root):
    with patch.object(subprocess, 'run', side_effect=AssertionError('native call')), \
         patch.object(subprocess, 'Popen', side_effect=AssertionError('native call')):
        return pkg.package(root, 'plan.json', 'approval.json', allow_synthetic_fixture=True)


def test_opaque_roundtrip_hydrates_every_file_and_publishes_receipt_last(tmp_path):
    from aav_benchmark.replay_workspace import hydrate_verified_archives
    plan = fixture(tmp_path); approve(tmp_path, plan); expected = pkg.inventory(tmp_path, plan)
    calls = []; original = pkg._publish_file
    def observed(source, *args):
        calls.append(source.name); return original(source, *args)
    with patch.object(pkg, '_publish_file', observed):
        receipt = run(tmp_path)
    assert calls[-1] == 'packaging_receipt.json'
    assert receipt['source_snapshot'] == expected and receipt['historical_first_job']['exact_original_paths_verified'] == 2
    assert receipt['accuracy_computed'] is receipt['numeric_pose_or_energy_payloads_interpreted'] is False
    hydrated = tmp_path / 'hydrated'; hydrated.mkdir(); entries = []
    for row in receipt['archives']:
        for kind in ['archive', 'manifest']:
            shutil.copyfile(tmp_path / plan['destination'] / row[kind]['path'], hydrated / row[kind]['path'])
        entries.append({k: row[k] for k in ['archive', 'manifest', 'target_prefix']})
    result = hydrate_verified_archives(hydrated, entries)
    assert result['verified_archive_members'] == len(expected['included'])
    for row in expected['included']:
        assert (hydrated / row['path']).read_bytes() == (tmp_path / row['path']).read_bytes()
    assert not (hydrated / 'tools').exists() and not (hydrated / 'data/processed').exists()


def test_exact_raw_and_entry_boundaries_and_single_file_refusal():
    row = lambda name, size: {'path': 'results/structural_pose_v3/' + name, 'bytes': size, 'sha256': '0' * 64}
    chunks = pkg.partition([row('c', 1), row('b', 0), row('a', 70_000_000)], pkg.LIMITS)
    assert [[r['bytes'] for r in c] for c in chunks] == [[70_000_000, 0], [1]]
    rows = [row(f'f{i:04d}', 0) for i in range(1001)]
    assert [len(c) for c in pkg.partition(list(reversed(rows)), pkg.LIMITS)] == [1000, 1]
    with pytest.raises(ValueError, match='single file'):
        pkg.partition([row('large', 70_000_001)], pkg.LIMITS)


def test_determinism_ignores_source_times_modes_and_traversal_order(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan); first = run(tmp_path)
    original = {pin['path']: (tmp_path / 'packaged' / pin['path']).read_bytes()
                for row in first['archives'] for pin in [row['archive'], row['manifest']]}
    for row in first['source_snapshot']['included']:
        path = tmp_path / row['path']; os.utime(path, (123456789, 123456789)); path.chmod(0o600)
    plan['destination'] = 'packaged_again'; approve(tmp_path, plan)
    walk = pkg.os.walk
    def reversed_walk(*args, **kwargs):
        for directory, dirs, files in walk(*args, **kwargs):
            yield directory, list(reversed(dirs)), list(reversed(files))
    with patch.object(pkg.os, 'walk', reversed_walk):
        second = run(tmp_path)
    assert second['source_snapshot_sha256'] == first['source_snapshot_sha256']
    assert {name: (tmp_path / 'packaged_again' / name).read_bytes() for name in original} == original


def test_finite_tree_and_cache_exclusions_hash_unknown_files_and_absent_future_tree(tmp_path):
    plan = fixture(tmp_path); tree = pkg.SOURCE_ROOTS['main'] + '/global_reports'
    cache = pkg.SOURCE_ROOTS['main'] + '/__pycache__/worker.pyc'
    put(tmp_path, tree + '/unknown/nested-name.json', b'opaque derived payload')
    put(tmp_path, cache, b'opaque cache bytes')
    plan['exclusions'] = [{'path': tree, 'scope': 'tree', 'kind': 'derived', 'reason': 'Separate derived publication'},
                          {'path': cache, 'kind': 'cache', 'reason': 'Interpreter cache'},
                          {'path': pkg.SOURCE_ROOTS['main'] + '/future-derived', 'scope': 'tree', 'kind': 'derived', 'reason': 'Future derived subtree'}]
    pkg.validate_plan(plan, allow_synthetic_fixture=True); approve(tmp_path, plan); receipt = run(tmp_path)
    excluded = receipt['source_snapshot']['excluded']
    assert {r['path'] for r in excluded} == {cache, tree + '/unknown/nested-name.json'}
    assert all(r['sha256'] == pkg.digest((tmp_path / r['path']).read_bytes()) for r in excluded)
    assert receipt['exclusion_rules'] == plan['exclusions']
    archived = {r['path'] for row in receipt['archives'] for r in json.loads((tmp_path / 'packaged' / row['manifest']['path']).read_text())['files']}
    assert not archived.intersection(r['path'] for r in excluded)


@pytest.mark.parametrize('change', ['approval', 'packager', 'plan_pin', 'writers', 'pending', 'accuracy', 'numeric', 'count', 'bool_count', 'unequal', 'snapshot'])
def test_root_approval_and_terminal_metadata_are_required(tmp_path, change):
    plan = fixture(tmp_path); ready_changes = {}; approval_changes = {}
    if change == 'approval': approval_changes['packaging_authorized'] = False
    elif change == 'packager': approval_changes['packager_sha256'] = '0' * 64
    elif change == 'plan_pin': approval_changes['packaging_plan'] = {'path': 'plan.json', 'bytes': 0, 'sha256': '0' * 64}
    elif change == 'writers': ready_changes['writers_stopped'] = False
    elif change == 'pending': ready_changes['pending'] = ['unfinished']
    elif change == 'accuracy': ready_changes['accuracy_computed'] = True
    elif change == 'numeric': ready_changes['numeric_pose_or_energy_payloads_interpreted'] = True
    elif change == 'snapshot': ready_changes['source_snapshot_sha256'] = '0' * 64
    else:
        counts = dict(pkg.READINESS_COUNTS, retained=2, recorded=2)
        if change == 'count': counts['main'] = 215
        elif change == 'bool_count': counts['retained'] = counts['recorded'] = True
        else: counts['recorded'] = 1
        ready_changes['counts'] = counts
    approve(tmp_path, plan, readiness_changes=ready_changes, approval_changes=approval_changes)
    with pytest.raises(ValueError):
        run(tmp_path)
    assert not (tmp_path / 'packaged').exists()


def test_unapproved_call_stops_before_any_source_inventory(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan, approval_changes={'status': 'proposal'})
    with patch.object(pkg, 'inventory', side_effect=AssertionError('source bytes opened')):
        with pytest.raises(ValueError, match='approval'):
            run(tmp_path)


def test_real_contract_and_cli_never_accept_synthetic_history(tmp_path):
    plan = fixture(tmp_path)
    with pytest.raises(ValueError, match='Real frozen'):
        pkg.validate_plan(plan)
    real = deepcopy(plan); real['profile'] = 'real_frozen_structural_outputs'; real['limits'] = dict(pkg.LIMITS)
    with pytest.raises(ValueError, match='All240'):
        pkg.validate_plan(real)
    real['historical_first_job'] = deepcopy(pkg.HISTORICAL)
    assert pkg.validate_plan(real) == real
    approve(tmp_path, plan)
    with pytest.raises(ValueError, match='Real frozen'):
        pkg.main(['--root', str(tmp_path), '--plan', 'plan.json', '--approval', 'approval.json'])


def test_existing_destination_and_coordinator_lock_are_preserved(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan)
    target = tmp_path / 'packaged'; target.mkdir(); (target / 'foreign').write_bytes(b'keep')
    with pytest.raises(FileExistsError): run(tmp_path)
    assert (target / 'foreign').read_bytes() == b'keep'
    shutil.rmtree(target)
    lock = tmp_path / pkg.MUST_BE_ABSENT[0]; lock.write_bytes(b'foreign coordinator')
    with pytest.raises(ValueError, match='lock must remain absent'): run(tmp_path)
    assert lock.read_bytes() == b'foreign coordinator' and not target.exists()


@pytest.mark.parametrize('kind', ['file_link', 'directory_link', 'root_link', 'fifo', 'installed', 'protected', 'worksharing', 'cache'])
def test_nonordinary_installed_and_protected_sources_are_rejected(tmp_path, kind):
    plan = fixture(tmp_path); base = tmp_path / pkg.SOURCE_ROOTS['main']
    outside = tmp_path / 'outside'; outside.mkdir(); (outside / 'opaque').write_bytes(b'outside')
    if kind == 'file_link': (base / 'linked').symlink_to(outside / 'opaque')
    elif kind == 'directory_link': (base / 'linked').symlink_to(outside, target_is_directory=True)
    elif kind == 'root_link':
        saved = tmp_path / 'saved'; base.rename(saved); base.symlink_to(saved, target_is_directory=True)
    elif kind == 'fifo': os.mkfifo(base / 'pipe')
    elif kind == 'installed': (base / 'site-packages').mkdir()
    elif kind == 'protected': (base / 'supplemental_arm').mkdir()
    elif kind == 'worksharing': put(tmp_path, pkg.SOURCE_ROOTS['main'] + '/worksharing_proposal.json', b'{}')
    else: put(tmp_path, pkg.SOURCE_ROOTS['main'] + '/__pycache__/source.pyc', b'cache')
    with pytest.raises(ValueError): pkg.inventory(tmp_path, plan)
    assert (outside / 'opaque').read_bytes() == b'outside'


@pytest.mark.parametrize('magic', [b'\x7fELF', b'MZ', b'\xcf\xfa\xed\xfe', b'\xca\xfe\xba\xbe', b'\x00asm'])
def test_executable_magic_cannot_be_hidden_as_derived_data(tmp_path, magic):
    plan = fixture(tmp_path); name = pkg.SOURCE_ROOTS['main'] + '/opaque.bin'
    put(tmp_path, name, magic + b'synthetic executable header')
    plan['exclusions'] = [{'path': name, 'kind': 'derived', 'reason': 'Must not waive executable prohibition'}]
    with pytest.raises(ValueError, match='executable'): pkg.inventory(tmp_path, plan)


@pytest.mark.parametrize('name', ['../outside', '/absolute', 'a//b', 'a/./b', 'a\\b', 'a:b', 'a/CON.txt', 'a/trailing.', 'a/\x00b'])
def test_nonportable_names_are_rejected(name):
    with pytest.raises(ValueError): pkg.portable_path(name)


def test_case_and_parent_collisions_are_platform_independent():
    for names in [[('a/File', 'file'), ('a/file', 'file')], [('a/DIR', 'file'), ('a/dir/child', 'file')],
                  [('a/empty', 'directory'), ('a/EMPTY', 'directory')]]:
        with pytest.raises(ValueError, match='collision'): pkg._namespace(names)


@pytest.mark.parametrize('when,operation', [('first', 'change'), ('last', 'change'), ('last', 'add'), ('last', 'delete'), ('last', 'excluded_add')])
def test_complete_snapshots_detect_all_source_changes_without_receipt(tmp_path, when, operation):
    plan = fixture(tmp_path); tree = pkg.SOURCE_ROOTS['main'] + '/global_reports'
    plan['exclusions'] = [{'path': tree, 'scope': 'tree', 'kind': 'derived', 'reason': 'Derived subtree'}]
    approve(tmp_path, plan); original = pkg._write_chunk
    def changed(root, stage, rows, number, total, *args):
        value = original(root, stage, rows, number, total, *args)
        if number == (1 if when == 'first' else total):
            path = root / pkg.SOURCE_ROOTS['main'] / 'history/partial.bin'
            if operation == 'change': path.write_bytes(b'changed source bytes')
            elif operation == 'delete': path.unlink()
            elif operation == 'add': put(root, pkg.SOURCE_ROOTS['assembly'] + '/new-file', b'new ordinary file')
            else: put(root, tree + '/new-derived.json', b'new opaque derived file')
        return value
    with patch.object(pkg, '_write_chunk', changed), pytest.raises(ValueError, match='changed|Source changed'):
        run(tmp_path)
    assert not (tmp_path / 'packaged/packaging_receipt.json').exists()
    assert not (tmp_path / 'packaged').exists()
    assert not list(tmp_path.glob('.structural-package-*'))


def test_during_read_hash_check_is_not_waived_by_restored_later_snapshot(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan); original = pkg._read; attempted = []
    source = pkg.SOURCE_ROOTS['main'] + '/history/partial.bin'
    def changed(root, name, **kwargs):
        raw = original(root, name, **kwargs)
        if name == source and kwargs.get('reject_executable'):
            attempted.append(name); return raw + b'not the snapshotted bytes'
        return raw
    with patch.object(pkg, '_read', changed), pytest.raises(ValueError, match='Source changed during'):
        run(tmp_path)
    assert attempted and not (tmp_path / 'packaged').exists()


def test_staged_archive_tampering_and_gzip_limit_abort_publication(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan); original = pkg._write_chunk
    def corrupt(root, stage, rows, number, total, *args):
        value = original(root, stage, rows, number, total, *args)
        (stage / value['archive']['path']).write_bytes(b'corrupted staged archive')
        return value
    with patch.object(pkg, '_write_chunk', corrupt), pytest.raises(ValueError, match='Pinned bytes changed'):
        run(tmp_path)
    assert not (tmp_path / 'packaged').exists()
    plan['limits']['max_gzip_bytes'] = 1; approve(tmp_path, plan)
    with pytest.raises(ValueError, match='exclusive size limit'): run(tmp_path)
    assert not (tmp_path / 'packaged').exists()


def test_cleanup_preserves_a_foreign_replacement(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan); original = pkg._publish_file; foreign = []
    def replaced(source, directory_fd, created, expected):
        original(source, directory_fd, created, expected)
        os.unlink(source.name, dir_fd=directory_fd)
        fd = os.open(source.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600, dir_fd=directory_fd)
        os.write(fd, b'foreign replacement'); os.close(fd); foreign.append(source.name)
        raise RuntimeError('synthetic publication interruption')
    with patch.object(pkg, '_publish_file', replaced), pytest.raises(RuntimeError): run(tmp_path)
    assert len(foreign) == 1
    assert (tmp_path / 'packaged' / foreign[0]).read_bytes() == b'foreign replacement'
    assert not (tmp_path / 'packaged/packaging_receipt.json').exists()


def test_normally_returning_publication_replacement_cannot_earn_receipt(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan); original = pkg._publish_file; foreign = []
    def replaced(source, directory_fd, created, expected):
        original(source, directory_fd, created, expected)
        if source.name.endswith('.tar.gz') and not foreign:
            os.unlink(source.name, dir_fd=directory_fd)
            fd = os.open(source.name, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600, dir_fd=directory_fd)
            os.write(fd, b'foreign replacement returned normally'); os.close(fd)
            foreign.append(source.name)
    with patch.object(pkg, '_publish_file', replaced), pytest.raises(ValueError, match='Published directory entry changed'):
        run(tmp_path)
    assert len(foreign) == 1
    assert (tmp_path / 'packaged' / foreign[0]).read_bytes() == b'foreign replacement returned normally'
    assert not (tmp_path / 'packaged/packaging_receipt.json').exists()


@pytest.mark.parametrize('change', ['contract', 'approval', 'lock'])
def test_late_contract_approval_or_stopped_writer_change_is_fatal(tmp_path, change):
    plan = fixture(tmp_path); approve(tmp_path, plan); original = pkg._write_chunk
    def changed(root, stage, rows, number, total, *args):
        value = original(root, stage, rows, number, total, *args)
        if number == total:
            if change == 'contract': (root / plan['source_contracts'][0]['path']).write_bytes(b'changed external source protocol')
            elif change == 'approval':
                approval = json.loads((root / 'approval.json').read_text()); approval['packaging_authorized'] = False
                write_json(root, 'approval.json', approval)
            else: (root / pkg.MUST_BE_ABSENT[0]).write_bytes(b'new coordinator lock')
        return value
    with patch.object(pkg, '_write_chunk', changed), pytest.raises(ValueError): run(tmp_path)
    assert not (tmp_path / 'packaged').exists()
    if change == 'lock': assert (tmp_path / pkg.MUST_BE_ABSENT[0]).read_bytes() == b'new coordinator lock'


@pytest.mark.parametrize('change', ['live_path', 'excluded', 'archive_member'])
def test_historical_original_paths_and_archive_members_remain_exact(tmp_path, change):
    plan = fixture(tmp_path); history = plan['historical_first_job']
    original = json.loads((tmp_path / history['manifest']['path']).read_text()); first = original['files'][0]
    if change == 'live_path': (tmp_path / first['path']).write_bytes(b'changed live historical bytes')
    elif change == 'excluded': plan['exclusions'] = [{'path': first['path'], 'kind': 'derived', 'reason': 'Not permitted for historical bytes'}]
    else:
        members = {row['path'][len(history['target_prefix']) + 1:]: (tmp_path / row['path']).read_bytes() for row in original['files']}
        history['archive'] = put(tmp_path, history['archive']['path'], historical_tar(members, corrupt=True))
    snapshot = pkg.inventory(tmp_path, plan)
    with pytest.raises(ValueError, match='Historical'):
        pkg.verify_historical(tmp_path, plan, snapshot)


@pytest.mark.parametrize('change', ['roots', 'destination', 'limit', 'wildcard', 'overlap', 'contracts', 'commit', 'must_absent'])
def test_plan_scope_and_limits_are_explicit(tmp_path, change):
    plan = fixture(tmp_path)
    if change == 'roots': plan['source_roots']['other'] = 'results/elsewhere'
    elif change == 'destination': plan['destination'] = pkg.SOURCE_ROOTS['main'] + '/archives'
    elif change == 'limit': plan['limits']['max_files'] = 1001
    elif change == 'wildcard': plan['exclusions'] = [{'path': pkg.SOURCE_ROOTS['main'] + '/*.json', 'kind': 'derived', 'reason': 'Invalid glob'}]
    elif change == 'overlap':
        path = pkg.SOURCE_ROOTS['main'] + '/derived'
        plan['exclusions'] = [{'path': path, 'scope': 'tree', 'kind': 'derived', 'reason': 'Tree'},
                              {'path': path + '/file', 'kind': 'derived', 'reason': 'Overlap'}]
    elif change == 'contracts': plan['source_contracts'][0]['path'] = 'tools/installed/source.json'
    elif change == 'commit': plan['source_protocol_commit'] = 'short'
    else: plan['must_be_absent'] = []
    with pytest.raises(ValueError): pkg.validate_plan(plan, allow_synthetic_fixture=True)


def test_outside_trees_are_not_opened_or_copied(tmp_path):
    plan = fixture(tmp_path)
    outside = ['tools/installed', 'data/processed/ignored', 'results/protected_supplemental/unrelated']
    for name in outside: put(tmp_path, name, b'\x7fELF never read or copied')
    approve(tmp_path, plan); receipt = run(tmp_path)
    assert all(not any(row['path'] == name or row['path'].startswith(name + '/') for name in outside)
               for row in receipt['source_snapshot']['included'])


def test_destination_claim_cannot_rebind_to_a_foreign_directory(tmp_path):
    plan = fixture(tmp_path); approve(tmp_path, plan); original = pkg.os.open
    destination = tmp_path / 'packaged'; replaced = []
    def opening(path, flags, *args, **kwargs):
        if Path(path) == destination and not replaced:
            destination.rename(tmp_path / 'retained-original-claim')
            destination.mkdir(); (destination / 'foreign').write_bytes(b'foreign directory')
            replaced.append(True)
        return original(path, flags, *args, **kwargs)
    with patch.object(pkg.os, 'open', opening), pytest.raises(ValueError, match='before opening'):
        run(tmp_path)
    assert (destination / 'foreign').read_bytes() == b'foreign directory'
    assert not (destination / 'packaging_receipt.json').exists()


@pytest.mark.parametrize('kind', ['root', 'ancestor'])
@pytest.mark.parametrize('entrypoint', ['inventory', 'package'])
def test_supplied_source_root_aliases_are_rejected_before_resolving(tmp_path, kind, entrypoint):
    container = tmp_path / 'ordinary'; container.mkdir()
    root = container / 'source'; root.mkdir(); plan = fixture(root); approve(root, plan)
    alias = tmp_path / 'alias'
    if kind == 'root':
        alias.symlink_to(root, target_is_directory=True); supplied = alias
    else:
        alias.symlink_to(container, target_is_directory=True); supplied = alias / 'source'
    with patch.object(pkg, '_authorization', side_effect=AssertionError('aliased root reached authorization')):
        with pytest.raises(ValueError, match='Supplied source root traverses a symlink'):
            if entrypoint == 'inventory': pkg.inventory(supplied, plan)
            else: pkg.package(supplied, 'plan.json', 'approval.json', allow_synthetic_fixture=True)
    assert not (root / 'packaged').exists()
