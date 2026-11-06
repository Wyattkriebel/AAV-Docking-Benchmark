"""Standard-library-only, explicit-commit exports and verified archive hydration."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import tarfile
import tempfile
import unicodedata


RECEIPT = ".replay/provenance.json"
CHUNK = 1024 * 1024


def digest_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json(raw):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Duplicate JSON object key")
            result[key] = value
        return result
    return json.loads(raw, object_pairs_hook=unique)


def relative_path(name, *, allow_empty=False):
    if allow_empty and name == "":
        return ""
    if (not isinstance(name, str) or not name or "\\" in name or ":" in name
            or unicodedata.normalize("NFC", name) != name or any(ord(c) < 32 or ord(c) == 127 for c in name)):
        raise ValueError("Path must be a portable repository-relative name")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != name or name == ".":
        raise ValueError("Path must be canonical and repository-relative")
    if any(p.rstrip(". ") != p or p.casefold() == ".git" for p in path.parts):
        raise ValueError("Unsafe path component")
    reserved = {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)), *(f"lpt{i}" for i in range(1, 10))}
    if any(p.split(".")[0].casefold() in reserved for p in path.parts):
        raise ValueError("Reserved device path component")
    parts = [p.casefold() for p in path.parts]
    if parts[0] in {"tools", ".replay"} or parts[:2] == ["data", "processed"]:
        raise ValueError("Tools, processed trees and replay metadata cannot be copied")
    return name


def _sha(value):
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("Expected an explicit lowercase SHA256")
    return value


def _pin(record):
    if not isinstance(record, dict) or not {"path", "sha256"} <= set(record) or set(record) - {"path", "sha256", "bytes"}:
        raise ValueError("Archive and manifest pins require path, SHA256 and optional byte size")
    relative_path(record["path"])
    _sha(record["sha256"])
    if "bytes" in record and (type(record["bytes"]) is not int or record["bytes"] < 0):
        raise ValueError("Pinned byte size must be a nonnegative integer")


def validate_plan(plan):
    if (not isinstance(plan, dict) or set(plan) != {"schema_version", "commit", "archives"}
            or type(plan["schema_version"]) is not int or plan["schema_version"] != 1
            or not isinstance(plan["commit"], str) or re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", plan["commit"]) is None
            or not isinstance(plan["archives"], list)):
        raise ValueError("Use schema1 with an exact full Git commit and explicit archive list")
    seen = set()
    for entry in plan["archives"]:
        if not isinstance(entry, dict) or set(entry) != {"archive", "manifest", "target_prefix"}:
            raise ValueError("Each archive requires archive/manifest pins and an explicit target_prefix")
        _pin(entry["archive"]); _pin(entry["manifest"])
        relative_path(entry["target_prefix"], allow_empty=True)
        if entry["archive"]["path"] in seen:
            raise ValueError("Duplicate archive in packaging plan")
        seen.add(entry["archive"]["path"])


def _environment():
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(GIT_NO_REPLACE_OBJECTS="1", GIT_NO_LAZY_FETCH="1", GIT_CONFIG_NOSYSTEM="1", GIT_CONFIG_GLOBAL=os.devnull, GIT_OPTIONAL_LOCKS="0")
    return env


def _command(repository, *arguments):
    return ["git", "--no-replace-objects", "--no-pager", "-C", str(repository), *arguments]


def _git(repository, *arguments):
    result = subprocess.run(_command(repository, *arguments), env=_environment(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode:
        raise ValueError(f"Local Git {arguments[0]} failed: {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def _repository(path, commit):
    path = Path(path).resolve(strict=True)
    bare = _git(path, "rev-parse", "--is-bare-repository").strip() == b"true"
    root = path if bare else Path(_git(path, "rev-parse", "--show-toplevel").decode().strip()).resolve()
    config = subprocess.run(_command(root, "config", "--local", "--get-regexp", r"(^extensions\.partialclone$|\.promisor$)"),
                            env=_environment(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if config.returncode not in (0, 1):
        raise ValueError("Cannot inspect local Git object availability policy")
    objects = Path(_git(root, "rev-parse", "--git-path", "objects").decode().strip())
    objects = objects if objects.is_absolute() else root / objects
    if config.returncode == 0 or any((objects / "pack").glob("*.promisor")):
        raise ValueError("Partial/promisor repositories are unsupported; all Git objects must be local")
    actual = _git(root, "rev-parse", "--verify", commit + "^{commit}").decode().strip()
    if actual != commit:
        raise ValueError("Requested Git object is not the exact commit")
    tree = _git(root, "rev-parse", commit + "^{tree}").decode().strip()
    git_dir = Path(_git(root, "rev-parse", "--absolute-git-dir").decode().strip()).resolve()
    common = Path(_git(root, "rev-parse", "--git-common-dir").decode().strip())
    common = common.resolve() if common.is_absolute() else (root / common).resolve()
    return root, (git_dir, common), tree


def _paths(entries):
    """Reject case aliases and file/directory conflicts before creating files."""
    nodes = {}
    for name in entries:
        relative_path(name)
        parts = PurePosixPath(name).parts
        for i in range(1, len(parts) + 1):
            prefix = PurePosixPath(*parts[:i]).as_posix()
            kind = "file" if i == len(parts) else "directory"
            previous = nodes.get(prefix.casefold())
            if previous is not None:
                if previous != (prefix, kind) or kind == "file":
                    raise ValueError("Duplicate, case-colliding or file/directory path conflict")
            else:
                nodes[prefix.casefold()] = (prefix, kind)


def _tree(repository, commit):
    entries = []
    for row in _git(repository, "ls-tree", "-r", "-z", "--full-tree", commit).split(b"\0"):
        if not row:
            continue
        metadata, name = row.split(b"\t", 1)
        mode, kind, object_id = metadata.decode("ascii").split()
        if mode not in {"100644", "100755"} or kind != "blob":
            raise ValueError("Git export permits only regular files; links and submodules are unsupported")
        entries.append((relative_path(name.decode("utf-8")), mode, object_id))
    _paths(name for name, _, _ in entries)
    return sorted(entries)


def _export(repository, commit, destination):
    entries = _tree(repository, commit)
    records = []
    with tempfile.TemporaryFile() as errors:
        process = subprocess.Popen(_command(repository, "cat-file", "--batch"), env=_environment(),
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=errors)
        try:
            for name, mode, object_id in entries:
                process.stdin.write((object_id + "\n").encode()); process.stdin.flush()
                header = process.stdout.readline().decode("ascii").strip().split()
                if len(header) != 3 or header[:2] != [object_id, "blob"] or not header[2].isdigit():
                    raise ValueError("Git batch returned an unexpected or unavailable blob")
                size = int(header[2]); remaining = size
                content_sha = hashlib.sha256()
                object_sha = hashlib.new("sha1" if len(object_id) == 40 else "sha256")
                object_sha.update(f"blob {size}\0".encode())
                path = destination / name; path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as output:
                    while remaining:
                        chunk = process.stdout.read(min(CHUNK, remaining))
                        if not chunk:
                            raise ValueError("Git blob was truncated")
                        remaining -= len(chunk); content_sha.update(chunk); object_sha.update(chunk); output.write(chunk)
                if process.stdout.read(1) != b"\n" or object_sha.hexdigest() != object_id:
                    raise ValueError("Exported Git blob identity changed")
                path.chmod(0o755 if mode == "100755" else 0o644)
                records.append({"path": name, "git_blob": object_id, "git_mode": mode, "bytes": size, "sha256": content_sha.hexdigest()})
            process.stdin.close(); process.stdout.close()
            if process.wait():
                raise ValueError("Git blob export failed")
        finally:
            if process.poll() is None:
                process.terminate(); process.wait()
    return records


def _check_file(root, record):
    path = root / record["path"]
    if path.is_symlink() or not path.is_file() or digest_file(path) != record["sha256"]:
        raise ValueError("Pinned archive/manifest is missing or has the wrong SHA256")
    if "bytes" in record and path.stat().st_size != record["bytes"]:
        raise ValueError("Pinned archive/manifest has the wrong size")
    return path


def _hydrate(root, entry, existing):
    archive_path = _check_file(root, entry["archive"])
    manifest_path = _check_file(root, entry["manifest"])
    manifest = _json(manifest_path.read_bytes())
    if (not isinstance(manifest, dict) or not isinstance(manifest.get("files"), list)
            or manifest.get("archive_sha256") != entry["archive"]["sha256"]
            or ("archive" in manifest and manifest["archive"] != archive_path.name)):
        raise ValueError("Manifest does not bind the declared archive")
    expected = {}
    for row in manifest["files"]:
        if not isinstance(row, dict) or not {"path", "bytes", "sha256"} <= set(row):
            raise ValueError("Every manifest member requires path, byte size and SHA256")
        name = relative_path(row["path"])
        _sha(row["sha256"])
        if name in expected or type(row["bytes"]) is not int or row["bytes"] < 0:
            raise ValueError("Duplicate or invalid manifest member")
        target = (PurePosixPath(entry["target_prefix"]) / name).as_posix()
        relative_path(target)
        expected[name] = {"path": target, "bytes": row["bytes"], "sha256": row["sha256"]}
    _paths(r["path"] for r in expected.values())
    # Check aliases against the complete exported/previously hydrated namespace.
    for row in expected.values():
        if row["path"].casefold() in {p.casefold() for p in existing} and row["path"] not in existing:
            raise ValueError("Archive target aliases an existing file path")
    _paths(set(existing) | {r["path"] for r in expected.values()})
    seen, records = set(), []
    with tarfile.open(archive_path, mode="r:*") as archive:
        for member in archive:
            name = relative_path(member.name)
            if (member.type not in (tarfile.REGTYPE, tarfile.AREGTYPE) or member.linkname
                    or getattr(member, "sparse", None) is not None):
                raise ValueError("Only ordinary regular archive members are allowed")
            if name in seen or name not in expected:
                raise ValueError("Duplicate or extra archive member")
            seen.add(name); row = expected[name]; path = root / row["path"]
            if member.size != row["bytes"]:
                raise ValueError("Archive member size differs from its manifest")
            reused = path.exists()
            if reused and (path.is_symlink() or not path.is_file() or path.stat().st_size != row["bytes"] or digest_file(path) != row["sha256"]):
                raise ValueError("Archive conflicts with an existing exported or hydrated file")
            path.parent.mkdir(parents=True, exist_ok=True)
            digest = hashlib.sha256(); size = 0
            output = None if reused else path.open("xb")
            try:
                source = archive.extractfile(member)
                for chunk in iter(lambda: source.read(CHUNK), b""):
                    size += len(chunk); digest.update(chunk)
                    if output is not None:
                        output.write(chunk)
            finally:
                if output is not None:
                    output.close()
            if size != row["bytes"] or digest.hexdigest() != row["sha256"]:
                raise ValueError("Archive member SHA256 or size differs from its manifest")
            origin = existing.get(row["path"], "created")
            records.append({"archive_member": name, **row, "disposition": origin + "_identical" if reused else "created"})
            existing.setdefault(row["path"], "archive")
    if seen != set(expected):
        raise ValueError("Archive is missing declared members")
    return {**entry, "members": records}


def _publish(stage, destination, *, existing=None):
    """Exclusively claim the destination; never rename over an appearing directory."""
    if existing is None:
        destination.mkdir(mode=0o700)
    def owned(path, digest=None):
        record = path.stat()
        return path, digest, (record.st_dev, record.st_ino)
    created = [owned(destination)] if existing is None else []
    try:
        files = sorted(p for p in stage.rglob("*") if p.is_file())
        receipt = stage / RECEIPT
        files = [p for p in files if p != receipt]
        if receipt.is_file():
            files.append(receipt)
        if existing is not None:
            files = [p for p in files if p.relative_to(stage).as_posix() not in existing]
        for source in files:
            target = destination / source.relative_to(stage)
            missing = []
            parent = target.parent
            while parent != destination and not parent.exists():
                missing.append(parent); parent = parent.parent
            for directory in reversed(missing):
                directory.mkdir(mode=0o700); created.append(owned(directory))
            if target.parent.is_symlink() or not target.parent.is_dir():
                raise ValueError("Publication parent changed to a non-directory")
            digest = digest_file(source)
            with source.open("rb") as stream, target.open("xb") as output:
                created.append(owned(target, digest))
                shutil.copyfileobj(stream, output, CHUNK)
                os.fchmod(output.fileno(), stat.S_IMODE(source.stat().st_mode))
            if digest_file(target) != digest:
                raise ValueError("Published file differs from verified staged bytes")
    except BaseException:
        # Remove only unchanged files created here, never a concurrent user's file.
        for path, digest, identity in reversed(created):
            try:
                current = path.lstat()
                if (current.st_dev, current.st_ino) != identity:
                    continue
                if digest is None:
                    path.rmdir()
                elif not path.is_symlink() and path.is_file() and digest_file(path) == digest:
                    path.unlink()
            except OSError:
                pass
        raise


def _git_ancestor(path):
    return any((p / ".git").exists() or (p / ".git").is_symlink()
               or ((p / "HEAD").is_file() and (p / "config").is_file() and (p / "objects").is_dir())
               for p in (path, *path.parents))


def hydrate_verified_archives(root, entries):
    """Hydrate a caller-owned private, newly assembled unversioned TEMP tree.

    This API records existing-file identity, not a Git-commit provenance claim.
    The caller must exclusively own this root while it is being assembled.
    """
    validate_plan({"schema_version": 1, "commit": "0" * 40, "archives": entries})
    requested = Path(root)
    if requested.is_symlink() or not requested.is_dir():
        raise ValueError("Hydration requires an ordinary caller-owned temporary directory")
    root = requested.resolve()
    if _git_ancestor(root):
        raise ValueError("Hydration may not write inside a Git working tree")
    files = []
    for path in root.rglob("*"):
        if path.is_symlink() or (not path.is_dir() and not path.is_file()):
            raise ValueError("Temporary hydration roots may contain only ordinary files/directories")
        if path.is_file():
            relative_path(path.relative_to(root).as_posix())
            files.append(path)
    for protected in (root / "tools", root / "data/processed", root / ".replay"):
        if protected.exists():
            raise ValueError("Temporary hydration roots must exclude protected trees and replay metadata")
    _paths(path.relative_to(root).as_posix() for path in files)
    existing = {path.relative_to(root).as_posix(): "existing" for path in files}
    with tempfile.TemporaryDirectory(prefix=".replay-stage-", dir=root.parent) as temporary:
        stage = Path(temporary)
        for source in files:
            target = stage / source.relative_to(root); target.parent.mkdir(parents=True, exist_ok=True)
            with source.open("rb") as stream, target.open("xb") as output:
                shutil.copyfileobj(stream, output, CHUNK)
        archives = [_hydrate(stage, entry, existing) for entry in entries]
        # Recheck every reused target before publishing any new file.
        for archive in archives:
            for member in archive["members"]:
                if member["disposition"] == "existing_identical":
                    _check_file(root, member)
        _publish(stage, root, existing={path.relative_to(root).as_posix() for path in files})
    return {"schema_version": 1, "archives": archives,
            "verified_archive_members": sum(len(a["members"]) for a in archives),
            "created_members": sum(r["disposition"] == "created" for a in archives for r in a["members"])}


def restore_workspace(repository, plan_path, destination):
    """Create one fresh outside-repository workspace; invoke no exported programs."""
    raw = Path(plan_path).read_bytes(); plan = _json(raw); validate_plan(plan)
    repository, git_dirs, tree = _repository(repository, plan["commit"])
    requested = Path(destination).absolute()
    if requested.exists() or requested.is_symlink():
        raise FileExistsError("Replay destination must not exist")
    parent = requested.parent.resolve(strict=True)
    destination = parent / requested.name
    if (_git_ancestor(parent) or destination.is_relative_to(repository)
            or any(destination.is_relative_to(p) for p in git_dirs)):
        raise ValueError("Replay workspace must remain outside the live repository")
    with tempfile.TemporaryDirectory(prefix=".replay-stage-", dir=parent) as temporary:
        stage = Path(temporary)
        exported = _export(repository, plan["commit"], stage)
        existing = {row["path"]: "git" for row in exported}
        archives = [_hydrate(stage, entry, existing) for entry in plan["archives"]]
        receipt = {"schema_version": 1, "commit": plan["commit"], "git_tree": tree,
                   "plan_sha256": hashlib.sha256(raw).hexdigest(), "plan": plan,
                   "exported_files": exported, "archives": archives}
        path = stage / RECEIPT; path.parent.mkdir()
        path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        _publish(stage, destination)
    return {"workspace": str(destination), "receipt": RECEIPT, "commit": plan["commit"],
            "exported_files": len(exported), "verified_archive_members": sum(len(a["members"]) for a in archives),
            "created_members": sum(r["disposition"] == "created" for a in archives for r in a["members"])}
