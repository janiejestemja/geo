# Quickstart

---

```bash
sudo dnf install python3-devel
```

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
python -m build coriolis_tools -n -w
pip install coriolis_tools/dist/*.whl
```
