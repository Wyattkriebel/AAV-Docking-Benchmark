"""Local Git objects and fresh temporary trees; no native or metric execution."""
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile

import pytest

from aav_benchmark import replay_workspace as replay


ROOT = Path(__file__).resolve().parents[1]


def encoded(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def git(repository, *arguments, data=None):
    env = replay._environment()
    env.update(GIT_AUTHOR_NAME="Replay test", GIT_AUTHOR_EMAIL="replay@example.invalid",
               GIT_COMMITTER_NAME="Replay test", GIT_COMMITTER_EMAIL="replay@example.invalid",
               GIT_AUTHOR_DATE="2026-11-03T00:00:00Z", GIT_COMMITTER_DATE="2026-11-03T00:00:00Z")
    return subprocess.check_output(["git", "-C", str(repository), *arguments], input=data, env=env).strip()


def commit_files(repository, files, *, object_format="sha1"):
    """Create an isolated object database through plumbing; never stage/checkout."""
    repository.mkdir(); git(repository, "init", "--quiet", "--object-format=" + object_format)
    tree = {}
    for name, value in files.items():
        mode, content = value if isinstance(value, tuple) else ("100644", value)
        node = tree; parts = name.split("/")
        for component in parts[:-1]: node = node.setdefault(component, {})
        oid = b"1" * 40 if mode == "160000" else git(repository, "hash-object", "-w", "--stdin", data=content)
        node[parts[-1]] = (mode, oid)
    def write_tree(node):
        rows = []
        for name, value in sorted(node.items()):
            if isinstance(value, dict): mode, kind, oid = "040000", "tree", write_tree(value)
            else:
                mode, oid = value; kind = "commit" if mode == "160000" else "blob"
            rows.append(f"{mode} {kind} ".encode() + oid + b"\t" + name.encode() + b"\n")
        return git(repository, "mktree", "--missing", data=b"".join(rows))
    return git(repository, "commit-tree", write_tree(tree).decode(), "-m", "fixture").decode()


def archive(entries):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as bundle:
        for name, content, kind in entries:
            member = tarfile.TarInfo(name); member.type = kind; member.size = len(content)
            if kind in (tarfile.SYMTYPE, tarfile.LNKTYPE): member.linkname = "outside"
            bundle.addfile(member, io.BytesIO(content))
    raw = stream.getvalue()
    return raw, {"archive": "inputs.tar.gz", "archive_sha256": digest(raw), "files": [
        {"path": name, "bytes": len(content), "sha256": digest(content)} for name, content, _ in entries]}


def fixture(tmp_path, entries=None, *, change=None, extra_git=None, prefix=""):
    entries = entries if entries is not None else [("shared.txt", b"same", tarfile.REGTYPE), ("nested/input.bin", b"new", tarfile.REGTYPE)]
    raw, manifest = archive(entries)
    if change: change(manifest)
    manifest_raw = encoded(manifest)
    files = {"shared.txt": b"same", "bundle/inputs.tar.gz": raw, "bundle/manifest.json": manifest_raw,
             ".gitattributes": b"shared.txt export-ignore\n*.txt export-subst\n"}
    files.update(extra_git or {})
    repository = tmp_path / "repository"; commit = commit_files(repository, files)
    plan = {"schema_version": 1, "commit": commit, "archives": [{
        "archive": {"path": "bundle/inputs.tar.gz", "sha256": digest(raw), "bytes": len(raw)},
        "manifest": {"path": "bundle/manifest.json", "sha256": digest(manifest_raw), "bytes": len(manifest_raw)}, "target_prefix": prefix}]}
    path = tmp_path / "plan.json"; path.write_bytes(encoded(plan))
    return repository, path, tmp_path / "workspace", plan


def test_exact_commit_export_ignores_worktree_and_attributes_and_reuses_identical_files(tmp_path):
    repository, plan_path, destination, plan = fixture(tmp_path, extra_git={"stamp.txt": b"$Format:%H$\n"})
    (repository / "shared.txt").write_bytes(b"uncommitted contents")
    (repository / "untracked-secret").write_bytes(b"do not copy")
    result = replay.restore_workspace(repository, plan_path, destination)
    assert result["commit"] == plan["commit"] and result["verified_archive_members"] == 2 and result["created_members"] == 1
    assert (destination / "shared.txt").read_bytes() == b"same"
    assert (destination / "nested/input.bin").read_bytes() == b"new"
    assert (destination / "stamp.txt").read_bytes() == b"$Format:%H$\n"
    assert not (destination / "untracked-secret").exists() and not (destination / ".git").exists()
    assert (repository / "shared.txt").read_bytes() == b"uncommitted contents"
    receipt = json.loads((destination / replay.RECEIPT).read_text())
    assert receipt["plan_sha256"] == digest(plan_path.read_bytes())
    assert {r["disposition"] for r in receipt["archives"][0]["members"]} == {"git_identical", "created"}
    assert all(replay.digest_file(destination / r["path"]) == r["sha256"] for r in receipt["exported_files"])


@pytest.mark.parametrize("name", ["../outside", "/absolute", "a/../outside", "a//b", "a\\b", "C:/outside",
                                 ".git/config", ".replay/provenance.json", "tools/input", "data/processed/input", "NUL.txt"])
def test_unsafe_or_protected_member_names_never_create_a_destination(tmp_path, name):
    repository, plan, destination, _ = fixture(tmp_path, [(name, b"input", tarfile.REGTYPE)])
    with pytest.raises(ValueError): replay.restore_workspace(repository, plan, destination)
    assert not destination.exists() and not list(tmp_path.glob(".replay-stage-*"))


@pytest.mark.parametrize("kind", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.DIRTYPE, tarfile.FIFOTYPE])
def test_archive_links_and_nonregular_types_are_rejected(tmp_path, kind):
    repository, plan, destination, _ = fixture(tmp_path, [("member", b"", kind)])
    with pytest.raises(ValueError, match="regular"): replay.restore_workspace(repository, plan, destination)
    assert not destination.exists()


@pytest.mark.parametrize("mutation", ["duplicate_tar", "duplicate_manifest", "missing", "extra", "size", "hash", "archive_hash"])
def test_full_archive_inventory_and_integrity_are_required(tmp_path, mutation):
    entries = [("member", b"input", tarfile.REGTYPE)] * (2 if mutation == "duplicate_tar" else 1)
    def change(manifest):
        if mutation == "duplicate_tar": manifest["files"] = manifest["files"][:1]
        elif mutation == "duplicate_manifest": manifest["files"] *= 2
        elif mutation == "missing": manifest["files"].append({"path": "missing", "bytes": 0, "sha256": digest(b"")})
        elif mutation == "extra": manifest["files"][0]["path"] = "different"
        elif mutation == "size": manifest["files"][0]["bytes"] += 1
        elif mutation == "hash": manifest["files"][0]["sha256"] = "0" * 64
        else: manifest["archive_sha256"] = "0" * 64
    repository, plan, destination, _ = fixture(tmp_path, entries, change=change)
    with pytest.raises(ValueError): replay.restore_workspace(repository, plan, destination)
    assert not destination.exists() and not list(tmp_path.glob(".replay-stage-*"))


@pytest.mark.parametrize("entries,extra", [
    ([("shared.txt", b"different", tarfile.REGTYPE)], {}),
    ([("shared.txt/child", b"new", tarfile.REGTYPE)], {}),
    ([("dir/child", b"new", tarfile.REGTYPE)], {"Dir/other": b"old"}),
    ([("Shared.txt", b"same", tarfile.REGTYPE)], {}),
])
def test_conflicting_or_case_aliasing_git_files_are_not_overwritten(tmp_path, entries, extra):
    repository, plan, destination, _ = fixture(tmp_path, entries, extra_git=extra)
    with pytest.raises(ValueError): replay.restore_workspace(repository, plan, destination)
    assert not destination.exists()


@pytest.mark.parametrize("mode", ["120000", "160000"])
def test_git_links_and_submodules_are_rejected(tmp_path, mode):
    repository, plan, destination, _ = fixture(tmp_path, extra_git={"link": (mode, b"outside")})
    with pytest.raises(ValueError, match="regular files"): replay.restore_workspace(repository, plan, destination)
    assert not destination.exists()


def test_fresh_destination_and_outside_repository_are_required(tmp_path):
    repository, plan, destination, _ = fixture(tmp_path)
    destination.mkdir(); sentinel = destination / "sentinel"; sentinel.write_bytes(b"preserve")
    with pytest.raises(FileExistsError): replay.restore_workspace(repository, plan, destination)
    assert sentinel.read_bytes() == b"preserve"
    with pytest.raises(ValueError, match="outside"): replay.restore_workspace(repository, plan, repository / "new-workspace")
    with pytest.raises(ValueError, match="outside"): replay.restore_workspace(repository, plan, repository / ".git/new-workspace")


def test_destination_appearing_during_publication_is_preserved(tmp_path, monkeypatch):
    repository, plan, destination, _ = fixture(tmp_path)
    original = replay._publish
    def race(stage, target):
        target.mkdir(); (target / "concurrent").write_bytes(b"preserve")
        return original(stage, target)
    monkeypatch.setattr(replay, "_publish", race)
    with pytest.raises(FileExistsError): replay.restore_workspace(repository, plan, destination)
    assert (destination / "concurrent").read_bytes() == b"preserve"
    assert not list(tmp_path.glob(".replay-stage-*"))


def test_explicit_prefix_and_offline_only_object_database(tmp_path):
    repository, plan, destination, _ = fixture(tmp_path, [("input", b"new", tarfile.REGTYPE)], prefix="results/example")
    replay.restore_workspace(repository, plan, destination)
    receipt = json.loads((destination / replay.RECEIPT).read_bytes())
    assert receipt["archives"][0]["members"][0]["path"] == "results/example/input"
    git(repository, "config", "remote.origin.promisor", "true")
    with pytest.raises(ValueError, match="Partial/promisor"): replay.restore_workspace(repository, plan, tmp_path / "partial")


def test_public_hydration_api_preserves_existing_bytes_and_rejects_git_roots(tmp_path):
    repository, _, _, plan = fixture(tmp_path)
    root = tmp_path / "private-source-tree"; root.mkdir()
    # The caller assembles only explicit source files before hydration.
    for name in ("shared.txt", "bundle/inputs.tar.gz", "bundle/manifest.json"):
        raw = subprocess.check_output(["git", "-C", str(repository), "show", plan["commit"] + ":" + name], env=replay._environment())
        path = root / name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(raw)
    before = (root / "shared.txt").read_bytes()
    result = replay.hydrate_verified_archives(root, plan["archives"])
    assert result["verified_archive_members"] == 2 and result["created_members"] == 1
    assert (root / "shared.txt").read_bytes() == before
    assert result["archives"][0]["members"][0]["disposition"] == "existing_identical"
    assert (root / "nested/input.bin").read_bytes() == b"new"
    with pytest.raises(ValueError, match="Git working tree"): replay.hydrate_verified_archives(repository, plan["archives"])


@pytest.mark.parametrize("change", ["schema_bool", "schema_float", "symbolic_commit", "command", "prefix", "pin_size"])
def test_packaging_plan_is_small_explicit_and_cannot_encode_commands(change):
    plan = {"schema_version": 1, "commit": "1" * 40, "archives": [{
        "archive": {"path": "inputs.tar.gz", "sha256": "a" * 64},
        "manifest": {"path": "manifest.json", "sha256": "b" * 64}, "target_prefix": ""}]}
    if change == "schema_bool": plan["schema_version"] = True
    elif change == "schema_float": plan["schema_version"] = 1.0
    elif change == "symbolic_commit": plan["commit"] = "HEAD"
    elif change == "command": plan["commands"] = ["arbitrary program"]
    elif change == "prefix": plan["archives"][0]["target_prefix"] = "../outside"
    else: plan["archives"][0]["archive"]["bytes"] = True
    with pytest.raises(ValueError): replay.validate_plan(plan)
    with pytest.raises(ValueError, match="Duplicate JSON"): replay._json(b'{"files":[],"files":[]}')


@pytest.mark.parametrize("field,key,value", [("archive", "sha256", "0" * 64), ("manifest", "sha256", "0" * 64), ("archive", "bytes", 1)])
def test_wrong_source_pins_fail_before_destination_publication(tmp_path, field, key, value):
    repository, path, destination, plan = fixture(tmp_path)
    plan["archives"][0][field][key] = value; path.write_bytes(encoded(plan))
    with pytest.raises(ValueError): replay.restore_workspace(repository, path, destination)
    assert not destination.exists() and not list(tmp_path.glob(".replay-stage-*"))


@pytest.mark.parametrize("extra", [{"tools/library": b"x"}, {"data/processed/frame": b"x"}, {".replay/content": b"x"}, {"Dir/a": b"a", "dir/b": b"b"}])
def test_protected_and_directory_alias_paths_in_git_are_rejected(tmp_path, extra):
    repository, path, destination, _ = fixture(tmp_path, extra_git=extra)
    with pytest.raises(ValueError): replay.restore_workspace(repository, path, destination)
    assert not destination.exists()


def test_verified_overlap_between_archives_is_reused_with_both_provenance_records(tmp_path):
    files, entries = {}, []
    for number in (1, 2):
        raw, manifest = archive([("shared.bin", b"same", tarfile.REGTYPE), (f"only{number}", b"new", tarfile.REGTYPE)])
        archive_name = f"bundle/input{number}.tar.gz"; manifest_name = f"bundle/manifest{number}.json"
        manifest["archive"] = Path(archive_name).name; manifest_raw = encoded(manifest)
        files.update({archive_name: raw, manifest_name: manifest_raw})
        entries.append({"archive": {"path": archive_name, "sha256": digest(raw)},
                        "manifest": {"path": manifest_name, "sha256": digest(manifest_raw)}, "target_prefix": "results/test"})
    repository = tmp_path / "repository"; commit = commit_files(repository, files)
    path = tmp_path / "plan.json"; path.write_bytes(encoded({"schema_version": 1, "commit": commit, "archives": entries}))
    destination = tmp_path / "workspace"; result = replay.restore_workspace(repository, path, destination)
    assert result["verified_archive_members"] == 4 and result["created_members"] == 3
    receipt = json.loads((destination / replay.RECEIPT).read_bytes())
    assert receipt["archives"][1]["members"][0]["disposition"] == "archive_identical"


def test_sha256_git_object_repository_exports_exact_bytes(tmp_path):
    repository = tmp_path / "repository"
    commit = commit_files(repository, {"input.txt": b"exact"}, object_format="sha256")
    assert len(commit) == 64
    plan = tmp_path / "plan.json"
    plan.write_bytes(encoded({"schema_version": 1, "commit": commit, "archives": []}))
    destination = tmp_path / "workspace"
    result = replay.restore_workspace(repository, plan, destination)
    assert result["commit"] == commit and (destination / "input.txt").read_bytes() == b"exact"


def test_other_linked_or_unrelated_git_worktrees_cannot_be_destinations(tmp_path):
    repository, plan_path, _, plan = fixture(tmp_path)
    git(repository, "update-ref", "refs/heads/main", plan["commit"])
    git(repository, "symbolic-ref", "HEAD", "refs/heads/main")
    linked = tmp_path / "linked"
    git(repository, "worktree", "add", "--detach", str(linked), plan["commit"])
    unrelated = tmp_path / "unrelated"; unrelated.mkdir(); git(unrelated, "init", "--quiet")
    for source, destination in ((linked, repository / "new-replay"), (repository, linked / "new-replay"),
                                (repository, unrelated / "new-replay")):
        with pytest.raises(ValueError, match="outside"):
            replay.restore_workspace(source, plan_path, destination)
        assert not destination.exists()


def test_published_698_member_preparation_archive_is_hydrated_exactly(tmp_path):
    manifest_path = ROOT / "results/prepared_input_foundation/prepared_inputs_manifest.json"
    manifest_raw = manifest_path.read_bytes(); manifest = json.loads(manifest_raw)
    archive_path = manifest_path.parent / manifest["archive"]; raw = archive_path.read_bytes()
    files = {manifest_path.relative_to(ROOT).as_posix(): manifest_raw, archive_path.relative_to(ROOT).as_posix(): raw}
    files.update({r["path"]: (ROOT / r["path"]).read_bytes() for r in manifest["indexes"]})
    repository = tmp_path / "repository"; commit = commit_files(repository, files)
    plan = {"schema_version": 1, "commit": commit, "archives": [{
        "archive": {"path": archive_path.relative_to(ROOT).as_posix(), "sha256": digest(raw)},
        "manifest": {"path": manifest_path.relative_to(ROOT).as_posix(), "sha256": digest(manifest_raw)}, "target_prefix": ""}]}
    plan_path = tmp_path / "plan.json"; plan_path.write_bytes(encoded(plan))
    destination = tmp_path / "workspace"
    result = replay.restore_workspace(repository, plan_path, destination)
    assert result["verified_archive_members"] == 698 and result["created_members"] == 696
    assert all((destination / r["path"]).stat().st_size == r["bytes"] and replay.digest_file(destination / r["path"]) == r["sha256"] for r in manifest["files"])
    assert not (destination / "tools").exists() and not (destination / "data/processed").exists()
