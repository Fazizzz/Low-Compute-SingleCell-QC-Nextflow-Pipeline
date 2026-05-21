process SOFTWARE_VERSIONS {
    label 'process_low'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path 'software_versions.yml', emit: yml

    script:
    """
    python3 - <<'PY' > software_versions.yml
import importlib, json, shutil, subprocess
def py_ver(mod):
    try:
        m = importlib.import_module(mod)
        return getattr(m, "__version__", "unknown")
    except Exception:
        return "not_installed"
def cmd_ver(cmd, flag="--version"):
    if shutil.which(cmd) is None:
        return "not_installed"
    try:
        out = subprocess.run([cmd, flag], capture_output=True, text=True, timeout=10)
        return (out.stdout or out.stderr).strip().split("\\n")[0]
    except Exception as e:
        return f"error:{e}"
data = {
    "kallisto": cmd_ver("kallisto", "version"),
    "bustools": cmd_ver("bustools", "version"),
    "kb_python": py_ver("kb_python"),
    "fastp": cmd_ver("fastp"),
    "scrublet": py_ver("scrublet"),
    "celltypist": py_ver("celltypist"),
    "kneed": py_ver("kneed"),
    "scanpy": py_ver("scanpy"),
    "anndata": py_ver("anndata"),
    "scipy": py_ver("scipy"),
    "numpy": py_ver("numpy"),
    "pandas": py_ver("pandas"),
    "plotly": py_ver("plotly"),
    "jinja2": py_ver("jinja2"),
}
for k, v in data.items():
    print(f"{k}: {v}")
PY
    """

    stub:
    """
    cat > software_versions.yml <<EOF
    kallisto: 0.51.1
    bustools: 0.44.1
    fastp: 1.3.3
    EOF
    """
}
