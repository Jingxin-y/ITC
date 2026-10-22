"""Verify every released file against the bundle's SHA256SUMS."""
from pathlib import Path
import hashlib

root = Path(__file__).resolve().parent
count = 0
for line in (root / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
    expected, relative = line.split("  ", 1)
    target = (root / relative).resolve()
    if not target.is_relative_to(root) or not target.is_file():
        raise SystemExit("Missing/invalid path: " + relative)
    if hashlib.sha256(target.read_bytes()).hexdigest() != expected:
        raise SystemExit("Checksum mismatch: " + relative)
    count += 1
print(f"PASS: {count} files match SHA256SUMS")
