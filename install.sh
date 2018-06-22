virtualenv venv/
source venv/bin/activate
pip install -r requirements_venv.txt
deactivate
(cd simulator && make)
(cd extractSequencesFromGTF/gffread && make)
