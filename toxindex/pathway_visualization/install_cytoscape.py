#!/usr/bin/env python3
"""
install_cytoscape.py

Cross‑platform download + install for Cytoscape v3.10.3.
Supports:
  --installer-path  : use a local file instead of downloading
  --target-dir      : install into this directory (defaults to system path)
  --dry-run         : show steps without executing
"""

import os
import sys
import platform
import argparse
import urllib.request
import subprocess
import shutil
import tarfile

VERSION    = "3.10.3"
BASE_URL   = "https://cytoscape.org"
DEFAULTS   = {
    "Darwin":   { "target": "/Applications", "archive": f"cytoscape-unix-{VERSION}.tar.gz" },
    "Linux":    { "target": "/opt",          "archive": f"cytoscape-unix-{VERSION}.tar.gz" },
    "Windows":  { "target": None,            "archive": f"Cytoscape_{VERSION}_windows_64bit.exe" }
}

def parse_args():
    p = argparse.ArgumentParser(
        description="Download & install Cytoscape v" + VERSION
    )
    p.add_argument(
        "--installer-path",
        help="Path to a pre‑downloaded installer/archive"
    )
    p.add_argument(
        "--target-dir",
        help="Override default install directory"
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print actions without executing"
    )
    return p.parse_args()

def detect_os():
    sys_name = platform.system()
    if sys_name not in DEFAULTS:
        sys.exit(f"Unsupported OS: {sys_name}")
    return sys_name

def get_installer(args, sys_name):
    info    = DEFAULTS[sys_name]
    archive = info["archive"]
    if args.installer_path:
        if not os.path.isfile(args.installer_path):
            sys.exit(f"Installer not found at {args.installer_path}")
        print(f"[+] Using local installer {args.installer_path}")
        return args.installer_path
    url   = f"{BASE_URL}/{archive}"
    print(f"[+] Downloading {archive} from {url}")
    urllib.request.urlretrieve(url, archive)
    print(f"[+] Saved to ./{archive}")
    return archive

def run(cmd, dry):
    print(f"→ {cmd}")
    if not dry:
        subprocess.run(cmd, shell=True, check=True)

def install_windows(installer, dry):
    # NSIS silent install is /S
    run(f'"{installer}" /S', dry)

def install_unix(installer, target_root, dry):
    tmp_dir = os.path.join("/tmp", f"cytoscape-{VERSION}")
    run(f"mkdir -p {tmp_dir}", dry)
    run(f"tar -xzf {installer} -C /tmp", dry)

    dest = os.path.join(target_root, f"Cytoscape-{VERSION}")
    run(f"rm -rf {dest}", dry)
    run(f"mv {tmp_dir} {dest}", dry)

    # symlink for easy CLI launch
    link = "/usr/local/bin/cytoscape"
    script = os.path.join(dest, "Cytoscape.sh")
    run(f"rm -f {link}", dry)
    run(f"ln -s {script} {link}", dry)

def main():
    args    = parse_args()
    sys_name= detect_os()

    # Determine where to install
    target_root = args.target_dir \
        or DEFAULTS[sys_name]["target"]

    # On Unix, if we're going system‑wide and not dry‑run, require sudo
    if sys_name in ("Darwin","Linux") \
       and not args.dry_run \
       and target_root in ("/opt","/Applications") \
       and os.geteuid()!=0:
        sys.exit("⚠ Please run with sudo or specify --target-dir to a writable path")

    installer = get_installer(args, sys_name)

    if sys_name == "Windows":
        install_windows(installer, args.dry_run)
    else:
        install_unix(installer, target_root, args.dry_run)

    print("[+] Done.")

if __name__=="__main__":
    main()
