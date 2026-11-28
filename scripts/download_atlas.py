import requests
import concurrent.futures
import os
import time
import urllib3
urllib3.disable_warnings()

url = "https://zenodo.org/records/21135969/files/atlas_celltype_annotated_481k.h5ad?download=1"
outfile = "data/atlas_celltype_annotated_481k.h5ad"
n_threads = 16

# get total size
head = requests.head(url, verify=False, allow_redirects=True, timeout=60)
total_size = int(head.headers.get("Content-Length", 0))
if total_size == 0:
    r = requests.get(url, headers={"Range": "bytes=0-0"}, verify=False, timeout=60)
    total_size = int(r.headers.get("Content-Range", "").split("/")[-1])
print("total size: %.2f GB" % (total_size / 1e9), flush=True)

chunk = total_size // n_threads

def download_segment(i):
    start = i * chunk
    end = (start + chunk - 1) if i < n_threads - 1 else (total_size - 1)
    partfile = "%s.part%d" % (outfile, i)
    want = end - start + 1
    if os.path.exists(partfile) and os.path.getsize(partfile) == want:
        return (i, want)  # already complete (resume)
    headers = {"Range": "bytes=%d-%d" % (start, end)}
    got = 0
    for attempt in range(3):  # retry up to 3 times
        try:
            r = requests.get(url, headers=headers, verify=False, stream=True, timeout=120)
            with open(partfile, "wb") as f:
                for block in r.iter_content(chunk_size=1024 * 1024):
                    if block:
                        f.write(block)
            got = os.path.getsize(partfile)
            if got == want:
                break
        except Exception as e:
            print("  segment %d attempt %d error: %s" % (i, attempt, e), flush=True)
            time.sleep(2)
    return (i, got)

t0 = time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=n_threads) as ex:
    results = list(ex.map(download_segment, range(n_threads)))

# verify all segments complete
missing = [r for r in results if r[1] != (chunk if r[0] < n_threads - 1 else total_size - (n_threads - 1) * chunk)]
if missing:
    print("INCOMPLETE segments:", missing, flush=True)
else:
    # merge
    with open(outfile, "wb") as out:
        for i in range(n_threads):
            partfile = "%s.part%d" % (outfile, i)
            with open(partfile, "rb") as f:
                out.write(f.read())
            os.remove(partfile)
    dt = time.time() - t0
    print("DONE: %.2f GB in %.1f min = %.2f MB/s" % (total_size / 1e9, dt / 60, total_size / dt / 1e6), flush=True)
