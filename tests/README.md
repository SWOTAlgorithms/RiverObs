# RiverObs testing steps:

First run RiverObs on the supplied RDF file:
```
estimate_swot_river.py riverobs-test.rdf
```
Next run pytest on the supplied suite of tests
```
pytest -v test_riverobs.py
```
...profit!
