"""
Downloads one real TNBC Hi-C sample (GSM5098082, "TNBC_Tissue3", the
smallest of GSE167150's three TNBC tumor .hic files) without fetching
the full GSE167150_RAW.tar (6.2GB, most of it unrelated cell-line and
normal-tissue files): GEO's plain HTTP directory listing only exposes
the combined tar, not per-sample files, but the tar's member order
matches its own filelist.txt manifest, so this walks the tar's
sequential 512-byte USTAR headers with tiny HTTP Range requests (each
header alone gives the exact byte size needed to jump to the next
header, with no data downloaded in between) until it finds the target
member, then issues exactly one Range GET for that member's ~645MB of
data -- a public technique for pulling one file out of a remote tar
without downloading the rest of it, not anything project-specific.
"""
import sys
import requests

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167150/suppl/GSE167150_RAW.tar"
BLOCK = 512


def parse_octal(field):
    field = field.rstrip(b"\x00").strip()
    return int(field, 8) if field else 0


def get_range(url, start, end):
    r = requests.get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=60)
    r.raise_for_status()
    return r.content


def walk_tar(url, target_name, max_entries=200):
    offset = 0
    for _ in range(max_entries):
        header = get_range(url, offset, offset + BLOCK - 1)
        if len(header) < BLOCK or header[0:1] == b"\x00":
            return None
        name = header[0:100].rstrip(b"\x00").decode("utf-8", "replace")
        size = parse_octal(header[124:136])
        data_start = offset + BLOCK
        if target_name in name:
            return dict(name=name, data_start=data_start, size=size)
        padded = ((size + BLOCK - 1) // BLOCK) * BLOCK
        offset = data_start + padded
    return None


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "GSM5098082_TNBC_Tissue3"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "../data/GSE167150_HiC/GSM5098082_TNBC_Tissue3.hic"

    print(f"Locating '{target}' inside the remote tar via Range-request header walk...")
    info = walk_tar(URL, target)
    if info is None:
        print(f"Target '{target}' not found.")
        sys.exit(1)
    print(f"Found {info['name']}: {info['size']} bytes at offset {info['data_start']}")

    start, size = info["data_start"], info["size"]
    end = start + size - 1
    print(f"Downloading via one Range GET (bytes={start}-{end}) to {out_path} ...")
    with requests.get(URL, headers={"Range": f"bytes={start}-{end}"},
                       stream=True, timeout=120) as r:
        r.raise_for_status()
        downloaded = 0
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                downloaded += len(chunk)
    print(f"Done: saved {downloaded} bytes to {out_path}")
