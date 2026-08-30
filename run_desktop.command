#!/bin/zsh
set -euo pipefail

APP_DIR="${0:A:h}"
APP_PYTHON="/opt/anaconda3/envs/rnafold/bin/python"

if [[ ! -x "$APP_PYTHON" ]]; then
  print -u2 "找不到 RNAfold 环境：$APP_PYTHON"
  print -u2 "请先确认 /opt/anaconda3/envs/rnafold 已安装。"
  read -r "?按回车键关闭…"
  exit 1
fi

export PYTHONNOUSERSITE=1
cd "$APP_DIR"
exec "$APP_PYTHON" -m primerfold.desktop
