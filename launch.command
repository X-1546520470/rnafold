#!/bin/zsh
# PrimerFold 浏览器版启动器（新电脑请先双击 setup.command 完成安装）
set -euo pipefail

APP_DIR="${0:A:h}"
cd "$APP_DIR"

candidates=(
  /opt/anaconda3/envs/rnafold/bin/python
  "$HOME/miniforge3/envs/rnafold/bin/python"
  "$HOME/miniconda3/envs/rnafold/bin/python"
  "$HOME/anaconda3/envs/rnafold/bin/python"
  /opt/miniforge3/envs/rnafold/bin/python
  /opt/miniconda3/envs/rnafold/bin/python
)
if command -v conda >/dev/null 2>&1; then
  candidates+=("$(conda info --base 2>/dev/null)/envs/rnafold/bin/python")
fi

for PY in $candidates; do
  if [[ -x "$PY" ]]; then
    export PYTHONNOUSERSITE=1
    exec "$PY" -m primerfold
  fi
done

print -u2 "未找到 rnafold conda 环境（需要 ViennaRNA）。"
print -u2 "请先双击同目录下的 setup.command 完成首次安装。"
read -r "?按回车键关闭…"
exit 1
