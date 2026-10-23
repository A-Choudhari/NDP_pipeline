#!/usr/bin/env python3
"""
test_federation.py -- one lab file, four federation fields, three standard formats.

    python3 tests/test_federation.py        # or: python3 -m pytest tests/

What these pin, 2026-10-01:
  * everything lab-specific comes from lab.yaml, so another lab edits one file;
  * every record says who owns the data, who may get it, whether it may train
    an AI model, and where the raw files are -- and is `incomplete` until it does;
  * every record is also written as schema.org/DCAT JSON-LD and as an RO-Crate,
    with the folder as one DCAT catalogue.

Fixtures are synthetic; nothing reads real data or writes outside a temp dir.
"""

import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
# The lab the tests run as, whatever the kit's own lab.yaml says.
os.environ["MD_PIPELINE_LAB"] = os.path.join(HERE, "lab_fixture.yaml")
os.environ.pop("NDP_ID_PREFIX", None)

import yaml                                      # noqa: E402

import ndp_build                                 # noqa: E402
from core import layout, ndp, settings, standard  # noqa: E402
from test_ndp import frame, md                   # noqa: E402

PIPELINE = os.path.dirname(HERE)


def read(path):
    with open(path) as fh:
        return fh.read()


def as_lab(lab_text, code):
    """Run `code` in a fresh interpreter whose lab file is `lab_text`, because
    the lab is read once, when the modules are imported."""
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "lab.yaml")
        with open(p, "w") as fh:
            fh.write(textwrap.dedent(lab_text))
        env = {**os.environ, "MD_PIPELINE_LAB": p}
        env.pop("NDP_ID_PREFIX", None)
        out = subprocess.run([sys.executable, "-c", textwrap.dedent(code)], cwd=PIPELINE,
                             env=env, capture_output=True, text=True)
        if out.returncode:
            raise AssertionError(out.stderr)
        return json.loads(out.stdout)


class TestLabFile(unittest.TestCase):
    PROBE = """
        import json
        from core import catalog, lab, ndp, topics
        print(json.dumps({
            "prefix": catalog.ID_PREFIX, "codes": catalog.TOPIC_CODE,
            "topic": catalog.topic_of("/data/KRAS_G12D", {}),
            "id": catalog.suggest_identifier("oncology", "/data/KRAS_G12D", "d/a"),
            "lab": ndp.DEFAULTS["lab"], "keywords": ndp.DEFAULTS["keywords"],
            "owner": ndp.DEFAULTS["data_owner"], "areas": list(topics.TOPIC_PATTERNS),
            "catalog": lab.info()["catalog_name"]}))
        """

    def test_another_lab_edits_one_file(self):
        got = as_lab("""
            lab: {name: "Smith Lab, Elsewhere U", catalog_name: Smith catalogue}
            identifiers: {prefix: smith}
            record_defaults:
              keywords: [molecular dynamics, Smith Lab]
              data_owner: {name: Jo Smith}
            topics:
              oncology: {code: ONC, words: [kras, p53]}
            index_topics:
              Oncology: ['\\bkras\\b']
            """, self.PROBE)
        self.assertEqual(got["prefix"], "SMITH")
        self.assertEqual(got["codes"], {"oncology": "ONC"})
        self.assertEqual(got["topic"], "oncology")
        self.assertTrue(got["id"].startswith("SMITH-ONC-KRAS"), got["id"])
        self.assertEqual(got["lab"], "Smith Lab, Elsewhere U")
        self.assertEqual(got["keywords"], ["molecular dynamics", "Smith Lab"])
        self.assertEqual(got["owner"]["name"], "Jo Smith")
        self.assertEqual(got["areas"], ["Oncology"])
        self.assertEqual(got["catalog"], "Smith catalogue")

    def test_an_empty_lab_file_names_no_lab(self):
        """Nothing in the code may fall back to one lab's values."""
        got = as_lab("{}\n", self.PROBE)
        self.assertEqual(got["prefix"], "LAB")
        self.assertIsNone(got["lab"])
        self.assertEqual(got["codes"], {})
        self.assertIsNone(got["topic"])
        self.assertIsNone(got["owner"]["name"])
        self.assertEqual(got["areas"], [])

    def test_a_section_left_empty_or_a_single_value_is_not_a_crash(self):
        got = as_lab("""
            lab:
              # nothing yet
            identifiers:
            record_defaults:
              keywords: Smith Lab
            storage:
              data_roots: /data/smith
            topics:
              onc: {code: ONC, words: kras}
            """, """
            import json
            from core import catalog, lab, ndp
            print(json.dumps({"prefix": catalog.ID_PREFIX, "kw": ndp.DEFAULTS["keywords"],
                              "roots": lab.data_roots(),
                              "topic": catalog.topic_of("/x/KRAS", {})}))
            """)
        self.assertEqual(got, {"prefix": "LAB", "kw": ["Smith Lab"],
                               "roots": ["/data/smith"], "topic": "onc"})

    def test_the_code_names_no_lab(self):
        """No lab's own values are in the code: the lab the worked examples
        were built for (examples/lab.yaml) is named nowhere in the kit's
        Python or shell files, only in lab files."""
        kit = os.path.dirname(PIPELINE)
        p = os.path.join(kit, "examples", "lab.yaml")
        if not os.path.exists(p):
            self.skipTest("no examples/lab.yaml")
        doc = yaml.safe_load(read(p)) or {}
        lab_block = doc.get("lab") or {}
        owner = (doc.get("record_defaults") or {}).get("data_owner") or {}
        needles = {lab_block.get("name"), lab_block.get("short_name"),
                   owner.get("name"), owner.get("email"),
                   f'"{(doc.get("identifiers") or {}).get("prefix")}"'}
        for account, who in (doc.get("people") or {}).items():
            needles |= {f'"{account}"', (who or {}).get("contact")}
            needles |= set((who or {}).get("authors") or [])
        needles |= set((doc.get("storage") or {}).get("data_roots") or [])
        needles = {n for n in needles if n and n != '"None"'}
        for root, dirs, files in os.walk(kit):
            dirs[:] = [d for d in dirs if d not in ("examples", "__pycache__", "tests")]
            for name in files:
                if not name.endswith((".py", ".sh")):
                    continue
                text = read(os.path.join(root, name))
                for n in needles:
                    self.assertNotIn(n, text, f"{os.path.relpath(os.path.join(root, name), kit)} "
                                              f"still says {n!r}")

class TestFederationFields(unittest.TestCase):
    def setUp(self):
        self.rows = frame(md("d", "1R", "a", source_root="/store/proj", source_dir="sys/1R"),
                          md("d", "2R", "a", source_root="/store/proj", source_dir="sys/2R"))
        self.m, _ = ndp.measure(self.rows, lambda *_: 1.01325)

    def build(self, cur=None, spec=None, defaults=None, warn=None):
        return ndp.build_record("d/a", self.rows, self.m, cur or {}, spec or {},
                                defaults or ndp.DEFAULTS, warn)

    def test_unanswered_until_someone_decides(self):
        rec = self.build()
        missing = ndp.missing_fields(rec)
        self.assertIn("access_level", missing)
        self.assertIn("ai_training", missing)
        # optional parts of a filled block are never what holds a record back
        self.assertNotIn("data_owner.orcid", missing)
        self.assertNotIn("data_location.uri", missing)

    def test_the_lab_answers_for_every_record(self):
        rec = self.build()
        self.assertEqual(rec["data_owner"]["name"], "Pat Investigator")
        self.assertEqual(rec["data_owner"]["email"], "pi@example.edu")
        self.assertEqual(rec["data_location"]["storage"], "Test Lab cluster storage")

    def test_where_the_raw_files_are_is_measured(self):
        rec = self.build()
        self.assertEqual(rec["data_location"]["path"], "/store/proj/sys")
        self.assertEqual([d["location"] for d in rec["sampling"]["replicas_detail"]],
                         ["/store/proj/sys/1R", "/store/proj/sys/2R"])

    def test_a_url_and_a_storage_name_come_from_different_levels(self):
        rec = self.build(spec={"data_location": {"uri": "https://data.example.org/d"}})
        self.assertEqual(rec["data_location"],
                         {"storage": "Test Lab cluster storage", "path": "/store/proj/sys",
                          "uri": "https://data.example.org/d"})

    def test_the_owner_is_never_assembled_from_two_people(self):
        """A system's own PI must not inherit the lab PI's e-mail."""
        rec = self.build(cur={"data_owner": "Dana Other"})
        self.assertEqual(rec["data_owner"],
                         {"name": "Dana Other", "email": None, "orcid": None,
                          "affiliation": None})

    def test_answers_are_read_the_way_people_write_them(self):
        for said, want in ((False, "prohibited"), (True, "allowed"), ("No", "prohibited"),
                           ("by agreement", "by-agreement"), ("ALLOWED", "allowed")):
            rec = self.build(cur={"ai_training": said, "access_level": "Licensable"})
            self.assertEqual(rec["ai_training"], want, said)
            self.assertEqual(rec["access_level"], "licensable")

    def test_an_answer_outside_the_list_is_refused_out_loud(self):
        warn = []
        rec = self.build(cur={"access_level": "semi-public", "ai_training": "maybe"},
                         warn=warn)
        self.assertIsNone(rec["access_level"])
        self.assertIsNone(rec["ai_training"])
        self.assertEqual(sum("is not one of" in w for w in warn), 2)

    def test_the_most_specific_level_wins(self):
        defaults = ndp.merged_defaults({"access_level": "private", "ai_training": None})
        self.assertEqual(defaults["access_level"], "private")
        rec = self.build(cur={"access_level": "public"},
                         spec={"access_level": "licensable", "ai_training": "no"},
                         defaults=defaults)
        self.assertEqual(rec["access_level"], "public")
        self.assertEqual(rec["ai_training"], "prohibited")

    def test_an_empty_key_in_settings_does_not_erase_the_lab(self):
        d = ndp.merged_defaults({"data_owner": None, "lab": None, "keywords": []})
        self.assertEqual(d["data_owner"]["name"], "Pat Investigator")
        self.assertEqual(d["lab"], "Test Lab, Example University")
        self.assertEqual(d["keywords"], [])      # an empty list is an answer
        # a block whose every key is empty is empty too
        rec = self.build(defaults={**ndp.DEFAULTS, "data_owner": {"name": None}},
                         spec={"data_owner": {"name": "", "email": None}})
        self.assertIsNone(rec["data_owner"]["name"])

    def test_an_empty_block_in_a_system_hides_nothing_from_its_dataset(self):
        """The build adds `data_location:` / `data_owner:` blocks with empty
        keys to a system; those must not hide what `ndp: datasets:` says."""
        doc = {"datasets": {"d": {"data_location": {"storage": "Globus X",
                                                    "uri": "globus://abc/d"},
                                  "data_owner": {"name": "Dana PI"}}},
               "systems": [{"id": "d/a", "data_location": {"path": None},
                            "data_owner": {"name": None}}]}
        rec = self.build(cur=ndp.curation_for("d/a", doc))
        self.assertEqual(rec["data_location"]["storage"], "Globus X")
        self.assertEqual(rec["data_location"]["uri"], "globus://abc/d")
        self.assertEqual(rec["data_owner"],
                         {"name": "Dana PI", "email": None, "orcid": None, "affiliation": None})
        # a system's own URL joins its dataset's storage, key by key
        doc["systems"][0]["data_location"] = {"uri": "https://own.example/a"}
        rec = self.build(cur=ndp.curation_for("d/a", doc))
        self.assertEqual((rec["data_location"]["storage"], rec["data_location"]["uri"]),
                         ("Globus X", "https://own.example/a"))

    def test_an_empty_block_in_defaults_leaves_the_lab_value(self):
        d = ndp.merged_defaults({"data_owner": {"name": None, "email": None},
                                 "data_location": {"storage": None, "uri": None}})
        self.assertEqual(d["data_owner"]["name"], "Pat Investigator")
        self.assertEqual(d["data_location"]["storage"], "Test Lab cluster storage")

    def test_replicas_on_different_filesystems_have_no_common_location(self):
        self.assertIsNone(ndp.raw_location([{"location": "/a/x"}, {"location": "/b/y"}]))
        self.assertEqual(ndp.raw_location([{"location": "/a/x/1"}, {"location": "/a/x/2"}]),
                         "/a/x")

    def test_new_defaults_survive_odd_keywords(self):
        text = ndp.scaffold_defaults([])
        old = ndp.DEFAULTS["keywords"]
        try:
            ndp.DEFAULTS["keywords"] = [2026, "HIV: Env"]
            doc = yaml.safe_load("ndp:\n" + ndp.scaffold_defaults([]))
        finally:
            ndp.DEFAULTS["keywords"] = old
        self.assertEqual(doc["ndp"]["defaults"]["keywords"], ["2026", "HIV: Env"])
        self.assertTrue(text)

    def test_the_catalogue_row_carries_them(self):
        from core import catalog
        rec = self.build(cur={"access_level": "private", "ai_training": "no"})
        row = catalog.entry(rec)
        self.assertEqual((row["data_owner"], row["access_level"], row["ai_training"]),
                         ("Pat Investigator", "private", "prohibited"))


class TestStandardFormats(unittest.TestCase):
    def setUp(self):
        rows = frame(md("d", "1R", "a", source_root="/store/proj", source_dir="sys/1R"),
                     md("d", "2R", "a", source_root="/store/proj", source_dir="sys/2R"))
        m, _ = ndp.measure(rows, lambda *_: 1.01325)
        self.rec = ndp.build_record(
            "d/a", rows, m,
            {"title": "BG505 trimer", "identifier": "TESTLAB-HIV-TEST-001", "topic": "hiv",
             "authors": ["A Researcher", "B Colleague"], "access_level": "licensable",
             "ai_training": "by-agreement"},
            {"contact": "x@example.edu", "license": "Creative Commons Attribution 4.0 International",
             "data_location": {"uri": "https://data.example.org/proj"}},
            ndp.DEFAULTS)
        self.entry = {"file": "sys_2x500ns.json", "keywords": ["molecular dynamics"],
                      "first_written": "2026-09-30", "last_changed": "2026-10-01"}

    def test_one_node_both_vocabularies(self):
        doc = standard.dataset_doc(self.rec, self.entry, "https://h.example/records")
        self.assertEqual(doc["@type"], ["Dataset", "dcat:Dataset"])
        self.assertEqual(doc["@id"], "https://h.example/records/sys_2x500ns.json#dataset")
        self.assertEqual(doc["name"], doc["dct:title"])
        self.assertEqual(doc["license"], "https://spdx.org/licenses/CC-BY-4.0")
        self.assertEqual(doc["dct:license"], {"@id": "https://spdx.org/licenses/CC-BY-4.0"})
        self.assertTrue(doc["dct:accessRights"]["@id"].endswith("/RESTRICTED"))
        self.assertFalse(doc["isAccessibleForFree"])
        self.assertEqual(doc["accountablePerson"]["name"], "Pat Investigator")
        self.assertEqual(doc["dct:rightsHolder"]["foaf:name"], "Pat Investigator")
        self.assertEqual([c["name"] for c in doc["creator"]], ["A Researcher", "B Colleague"])
        self.assertEqual(doc["dct:issued"], {"@value": "2026-09-30", "@type": "xsd:date"})
        raw = [d for d in doc["distribution"] if d["name"] == "Raw trajectories"][0]
        self.assertEqual(raw["contentUrl"], "https://data.example.org/proj")
        self.assertIn("/store/proj/sys", raw["description"])
        self.assertEqual({d["@id"] for d in doc["distribution"]},
                         {d["@id"] for d in doc["dcat:distribution"]})

    def test_ai_training_is_an_odrl_policy(self):
        doc = standard.dataset_doc(self.rec, self.entry)
        pol = doc["odrl:hasPolicy"]
        self.assertEqual(pol["@type"], "odrl:Set")
        rule = pol["odrl:permission"][0]
        self.assertEqual(rule["odrl:constraint"]["odrl:rightOperand"], standard.AI_PURPOSE)
        self.assertEqual(rule["odrl:duty"]["odrl:action"]["@id"], "odrl:obtainConsent")
        self.assertIn("agreement", doc["usageInfo"]["text"])
        no = standard.dataset_doc({**self.rec, "ai_training": "prohibited"}, self.entry)
        self.assertIn("odrl:prohibition", no["odrl:hasPolicy"])
        unset = standard.dataset_doc({**self.rec, "ai_training": None}, self.entry)
        self.assertNotIn("odrl:hasPolicy", unset)

    def test_replicas_whose_ids_slug_alike_stay_two(self):
        rec = json.loads(json.dumps(self.rec))
        rec["sampling"]["replicas_detail"][0]["id"] = "r 1"
        rec["sampling"]["replicas_detail"][1]["id"] = "r/1"
        graph = standard.crate(rec, self.entry)["@graph"]
        runs = [n for n in graph if n["@id"].startswith("#run-")]
        self.assertEqual(len(runs), 2)

    def test_two_records_never_share_an_id(self):
        a = standard.dataset(self.rec, {**self.entry, "file": "a_3x1us.json"})
        b = standard.dataset(self.rec, {**self.entry, "file": "b_3x1us.json"})
        self.assertNotEqual(a["@id"], b["@id"])        # same identifier, two files

    def test_a_contact_that_is_not_one_address_is_no_mailto(self):
        rec = {**self.rec, "contact": "a@x.edu, b@y.edu"}
        doc = standard.dataset(rec, self.entry)
        self.assertEqual(doc["dcat:contactPoint"]["vcard:hasEmail"], "a@x.edu, b@y.edu")
        crate = standard.crate(rec, self.entry)
        self.assertFalse(any(n["@id"].startswith("mailto:a@x.edu,") for n in crate["@graph"]))

    def test_an_engine_line_splits_into_name_and_version(self):
        rec = {"setup": {"software": "NAMD 2.14 (Frontera, TACC); AMBER 22; OpenMM"}}
        self.assertEqual(standard.software_of(rec),
                         [("NAMD", "2.14", "Frontera, TACC"), ("AMBER", "22", None),
                          ("OpenMM", None, None)])

    def test_nothing_null_is_written(self):
        rec = {**self.rec, "related_publication": None, "contact": None}
        text = json.dumps(standard.dataset_doc(rec, self.entry))
        self.assertNotIn("null", text)
        self.assertNotIn("citation", text)

    def test_the_crate_is_a_process_run_crate(self):
        crate = standard.crate(self.rec, self.entry)
        self.assertEqual(crate["@context"], "https://w3id.org/ro/crate/1.1/context")
        by_id = {n["@id"]: n for n in crate["@graph"]}
        self.assertEqual(by_id["ro-crate-metadata.json"]["about"], {"@id": "./"})
        root = by_id["./"]
        for key in ("name", "description", "datePublished", "license"):
            self.assertIn(key, root)
        self.assertEqual(root["conformsTo"], {"@id": standard.PROCESS_RUN})
        self.assertEqual({p["@id"] for p in root["hasPart"]}, {"record.json", "dataset.jsonld"})
        runs = [n for n in crate["@graph"] if n["@type"] == "CreateAction"]
        self.assertEqual(len(runs), 3)          # two replicas, one metadata read
        replica = by_id["#run-1R"]
        engine = by_id[replica["instrument"]["@id"]]
        self.assertEqual((engine["name"], engine["softwareVersion"]), ("NAMD", "3.0.2"))
        raw = by_id[replica["result"]["@id"]]
        self.assertEqual(raw["@id"], "https://data.example.org/proj/1R")
        for n in crate["@graph"]:               # every reference resolves
            for v in n.values():
                for ref in (v if isinstance(v, list) else [v]):
                    if isinstance(ref, dict) and set(ref) == {"@id"} \
                            and not ref["@id"].startswith(("http", "mailto")):
                        self.assertIn(ref["@id"], by_id, f"{n['@id']} -> {ref}")


class TestBuild(unittest.TestCase):
    def build(self, tmp, *rows, args=()):
        reg = layout.mddb_registry(tmp)
        os.makedirs(reg, exist_ok=True)
        frame(*rows).to_parquet(os.path.join(reg, "mds.parquet"), index=False)
        if not settings.exists(tmp):
            settings.create(tmp, folder_section=f"folder:\n  root: {tmp}\n")
        return ndp_build.main([tmp, *args])

    def test_a_build_writes_every_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            results = layout.results(tmp)
            rec_name = [f for f in os.listdir(results)
                        if f.endswith(".json") and f != "catalog.json"][0]
            d = os.path.join(results, standard.DIR, rec_name[:-5])
            self.assertEqual(sorted(os.listdir(d)), sorted(standard.OURS))
            self.assertEqual(read(os.path.join(d, "record.json")),
                             read(os.path.join(results, rec_name)))
            cat = json.loads(read(os.path.join(results, standard.DIR, standard.CATALOG)))
            self.assertEqual(len(cat["dcat:dataset"]), 1)
            rec = json.loads(read(os.path.join(results, rec_name)))
            for key in ("data_owner", "access_level", "ai_training", "data_location"):
                self.assertIn(key, rec)
            self.assertEqual(rec["data_location"]["path"], os.path.join(tmp, "1R"))
            sheet = read(layout.registration_md(tmp))
            self.assertIn("| Access level | **missing** |", sheet)
            self.assertIn("| Data owner (PI) | Pat Investigator", sheet)

            # the new settings.yaml asks the questions, and still parses
            doc = yaml.safe_load(read(settings.path(tmp)))
            dflt = doc["ndp"]["defaults"]
            self.assertIn("access_level", dflt)
            # left empty, so the folder keeps following lab.yaml
            self.assertIsNone(dflt["data_owner"])
            self.assertIn("Pat Investigator", read(settings.path(tmp)))   # in the comment

            # an unchanged build changes no file
            before = {f: os.path.getmtime(os.path.join(d, f)) for f in os.listdir(d)}
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            self.assertEqual(before,
                             {f: os.path.getmtime(os.path.join(d, f)) for f in os.listdir(d)})

    def test_answering_in_settings_makes_the_record_ready_for_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            p = settings.path(tmp)
            text = read(p)
            text = text.replace("  access_level:", "  access_level: public", 1)
            text = text.replace("  ai_training:", "  ai_training: allowed", 1)
            with open(p, "w") as fh:
                fh.write(text)
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            idx = json.loads(read(layout.ndp_index(tmp)))["records"][0]
            self.assertNotIn("access_level", idx["missing"])
            self.assertNotIn("ai_training", idx["missing"])

    def test_every_missing_field_has_a_line_in_its_own_block(self):
        """The place to answer a missing field is always in the system's
        block, even for a field that went missing after the block was
        written -- and nothing a person wrote is touched."""
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            p = settings.path(tmp)
            block = yaml.safe_load(read(p))["ndp"]["systems"][0]
            self.assertIn("access_level", block)      # added by the build
            self.assertIn("construct", block)
            # a person's text, and a key that has gone from the block
            text = read(p).replace("    title:\n", "    title: My trimer  # mine\n")
            text = text.replace("    construct:\n", "", 1)
            with open(p, "w") as fh:
                fh.write(text)
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            after = read(p)
            self.assertIn("    title: My trimer  # mine\n", after)
            block = yaml.safe_load(after)["ndp"]["systems"][0]
            self.assertIn("construct", block)
            self.assertIsNone(block["construct"])
            # nothing left to add: a rebuild leaves the file as it is
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            self.assertEqual(read(p), after)
            # an answer typed on the added line reaches the record
            with open(p, "w") as fh:
                fh.write(after.replace("    access_level:", "    access_level: private", 1))
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            idx = json.loads(read(layout.ndp_index(tmp)))["records"][0]
            self.assertNotIn("access_level", idx["missing"])

    def test_lines_go_in_whatever_the_file_looks_like(self):
        """Windows line endings stay; a `|` text with blank lines inside a
        block is still that block; another system's odd block does not stop
        this one from getting its lines."""
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            p = settings.path(tmp)
            text = read(p).replace("    construct:\n", "", 1).replace(
                "    description:\n", "    description: |\n      First.\n\n      Second.\n", 1)
            with open(p, "w", newline="") as fh:
                fh.write(text.replace("\n", "\r\n"))
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            with open(p, newline="") as fh:
                raw = fh.read()
            self.assertNotIn("\n", raw.replace("\r\n", ""))      # every line still CRLF
            block = yaml.safe_load(raw)["ndp"]["systems"][0]
            self.assertIn("construct", block)
            self.assertEqual(block["description"], "First.\n\nSecond.\n")

    def test_check_adds_no_lines(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            p = settings.path(tmp)
            text = read(p).replace("    construct:\n", "", 1)
            with open(p, "w") as fh:
                fh.write(text)
            self.build(tmp, md("d", "1R", "aaaa", source_root=tmp), args=["--check"])
            self.assertEqual(read(p), text)

    def test_the_sweep_removes_only_its_own_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            top = os.path.join(layout.results(tmp), standard.DIR)
            stale, mine = os.path.join(top, "gone_1x500ns"), os.path.join(top, "notes")
            for d, f in ((stale, standard.DATASET), (mine, "README.md")):
                os.makedirs(d)
                with open(os.path.join(d, f), "w") as fh:
                    fh.write("{}")
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            self.assertFalse(os.path.exists(stale))
            self.assertTrue(os.path.exists(os.path.join(mine, "README.md")))

    def test_check_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(self.build(tmp, md("d", "1R", "aaaa", source_root=tmp)), 0)
            top = os.path.join(layout.results(tmp), standard.DIR)
            import shutil
            shutil.rmtree(top)
            self.build(tmp, md("d", "1R", "aaaa", source_root=tmp), args=["--check"])
            self.assertFalse(os.path.exists(top))


if __name__ == "__main__":
    unittest.main(verbosity=2)
