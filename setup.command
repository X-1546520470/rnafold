#!/bin/zsh
# PrimerFold 首次安装脚本（新电脑上双击运行一次即可）
# 做三件事: 找到/提示安装 conda → 创建 rnafold 环境（含 ViennaRNA）→ 安装依赖
set -euo pipefail

APP_DIR="${0:A:h}"
cd "$APP_DIR"

print -u2 "PrimerFold 首次安装"
print -u2 "=================="

# ---- 1. 寻找 conda ----
CONDA=""
for c in \
  "$HOME/miniforge3/bin/conda" \
  "$HOME/miniconda3/bin/conda" \
  "$HOME/anaconda3/bin/conda" \
  /opt/miniforge3/bin/conda \
  /opt/miniconda3/bin/conda \
  /opt/anaconda3/bin/conda; do
  if [[ -x "$c" ]]; then CONDA="$c"; break; fi
done

if [[ -z "$CONDA" ]] && command -v conda >/dev/null 2>&1; then
  CONDA="$(command -v conda)"
fi

if [[ -z "$CONDA" ]]; then
  print -u2 ""
  print -u2 "未找到 conda。请先安装 Miniforge（免费、轻量）："
  print -u2 "  1. 打开 https://github.com/conda-forge/miniforge#download"
  print -u2 "  2. 下载 macOS 对应版本（Apple 芯片选 arm64）并安装"
  print -u2 "  3. 重新双击本脚本"
  read -r "?按回车键关闭…"
  exit 1
fi
print -u2 "找到 conda: $CONDA"

# ---- 2. 创建 rnafold 环境（含 ViennaRNA + Python + tkinter） ----
ENV_NAME="rnafold"
if "$CONDA" env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
  print -u2 "conda 环境 ${ENV_NAME} 已存在，跳过创建。"
else
  print -u2 "正在创建 conda 环境 ${ENV_NAME}（含 ViennaRNA，约需几分钟）…"
  "$CONDA" env create -n "$ENV_NAME" -f "$APP_DIR/environment.yml"
fi

ENV_PYTHON="$("$CONDA" info --base)/envs/${ENV_NAME}/bin/python"
if [[ ! -x "$ENV_PYTHON" ]]; then
  print -u2 "环境已创建但未找到 python: $ENV_PYTHON"
  read -r "?按回车键关闭…"
  exit 1
fi

# ---- 3. 安装 PrimerFold 依赖（可编辑模式） ----
print -u2 "安装 PrimerFold 依赖…"
PYTHONNOUSERSITE=1 "$ENV_PYTHON" -m pip install -e "$APP_DIR" -q

# ---- 4. 自检 ----
print -u2 "自检 ViennaRNA 命令…"
for tool in RNAfold RNAduplex RNAcofold RNAplot; do
  if ! "$ENV_PYTHON" -c "import shutil,sys; sys.exit(0 if shutil.which('$tool') or __import__('os').access('$(dirname "$ENV_PYTHON")/$tool', os.X_OK) else 1)" 2>/dev/null; then
    print -u2 "⚠ 未找到 $tool（ViennaRNA 可能未装全）"
  fi
done

print -u2 ""
print -u2 "✔ 安装完成！以后双击 run_desktop.command（桌面版）或"
print -u2 "  launch.command（浏览器版）即可使用。"
read -r "?按回车键关闭…"
