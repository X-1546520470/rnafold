"""Launch the PrimerFold GUI with ``python -m primerfold``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        prog="primerfold",
        description="PrimerFold：本地 DNA 引物发卡与二聚体分析 GUI（ViennaRNA 后端）。",
    )
    parser.add_argument(
        "--address",
        default="127.0.0.1",
        help="监听地址（默认 127.0.0.1，仅本机访问）。",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=8501,
        help="监听端口（默认 8501）。",
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="不自动打开浏览器。",
    )
    args = parser.parse_args()

    app_path = Path(__file__).with_name("app.py")
    try:
        from streamlit.web import cli as stcli
    except ImportError:
        print(
            "缺少 streamlit，请先安装依赖：pip install -e . 或 pip install streamlit pandas",
            file=sys.stderr,
        )
        return 1

    streamlit_args = [
        "streamlit",
        "run",
        str(app_path),
        f"--server.address={args.address}",
        f"--server.port={args.port}",
    ]
    if args.headless:
        streamlit_args.append("--server.headless=true")
    sys.argv = streamlit_args
    return stcli.main()


if __name__ == "__main__":
    raise SystemExit(main())
