"""从清华镜像抓取 cp313/win_amd64 wheel 到 tools/wheelhouse，供离线安装。

用法: .venv/Scripts/python.exe scripts/fetch_wheels.py
全部路径基于脚本自身所在目录推导，不含绝对路径。
"""
import re
import sys
import urllib.request
from urllib.parse import urljoin
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WHEELHOUSE = ROOT / "tools" / "wheelhouse"
WHEELHOUSE.mkdir(parents=True, exist_ok=True)

MIRROR = "https://pypi.tuna.tsinghua.edu.cn/simple/{}/"
UA = {"User-Agent": "curl/8.5.0"}

PKGS = [
    "numpy", "pandas", "scipy", "scikit-learn", "xgboost", "rdkit",
    "matplotlib", "joblib", "python-dateutil", "six", "pytz", "tzdata",
    "cycler", "kiwisolver", "pyparsing", "fonttools", "packaging",
    "pillow", "threadpoolctl", "contourpy",
]


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=90).read().decode("utf-8", "ignore")


def ver_key(fname):
    """从 wheel 文件名取版本号，转成可比较的数值元组。"""
    parts = fname.split("-")
    if len(parts) < 2:
        return ()
    nums = re.findall(r"\d+", parts[1])
    return tuple(int(n) for n in nums[:4])


def pick(html, base):
    links = re.findall(r'href="([^"#]+\.whl)(?:#[^"]*)?"', html)
    ok = []
    for l in links:
        f = l.rsplit("/", 1)[-1]
        if "cp313t" in f:
            continue
        if "win_amd64" in f and ("cp313" in f or "abi3" in f or "py3-none" in f or "py2.py3-none" in f):
            ok.append(urljoin(base, l))
        elif "none-any" in f:
            ok.append(urljoin(base, l))
    if not ok:
        return None
    ok.sort(key=lambda l: ver_key(l.rsplit("/", 1)[-1]))
    return ok[-1]


def main():
    pkgs = sys.argv[1:] or PKGS
    fail = []
    for pkg in pkgs:
        try:
            html = fetch(MIRROR.format(pkg))
            link = pick(html, MIRROR.format(pkg))
            if not link:
                print(f"[MISS] {pkg}: 无匹配 wheel")
                fail.append(pkg)
                continue
            fname = link.rsplit("/", 1)[-1]
            dest = WHEELHOUSE / fname
            if dest.exists() and dest.stat().st_size > 0:
                print(f"[SKIP] {fname}")
                continue
            data = urllib.request.urlopen(
                urllib.request.Request(link, headers=UA),
                timeout=180,
            ).read()
            dest.write_bytes(data)
            print(f"[OK]   {fname}  {len(data)/1e6:.1f} MB")
        except Exception as e:
            print(f"[ERR]  {pkg}: {type(e).__name__} {e}")
            fail.append(pkg)
    if fail:
        print("\n失败包:", ", ".join(fail))
        sys.exit(1)
    print("\n全部 wheel 就绪:", WHEELHOUSE)


if __name__ == "__main__":
    main()
