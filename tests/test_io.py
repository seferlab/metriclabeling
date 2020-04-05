import bz2
import gzip

from metriclabeling.io import IdMapper, read_biogrid, read_gaf, read_obo

OBO = """format-version: 1.2
default-namespace: gene_ontology

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0000001
name: a
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000002
name: b
namespace: biological_process
alt_id: GO:0000099
is_a: GO:0008150
relationship: part_of GO:0000001 ! a

[Term]
id: GO:0000003
name: c
namespace: biological_process
relationship: regulates GO:0000001

[Term]
id: GO:0000004
name: old
namespace: biological_process
is_obsolete: true

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component

[Typedef]
id: part_of
name: part of
"""


def test_obo(tmp_path):
    p = tmp_path / "go.obo.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(OBO)
    onto = read_obo(p, namespace="P")
    assert len(onto) == 4  # obsolete and CC dropped
    b = onto.index["GO:0000002"]
    assert sorted(onto.ids[x] for x in onto.parents[b]) == ["GO:0000001", "GO:0008150"]
    c = onto.index["GO:0000003"]
    assert onto.parents[c] == []  # 'regulates' not kept by default
    assert onto.resolve("GO:0000099") == b
    assert onto.up_distances(b)[onto.index["GO:0008150"]] == 1
    assert set(onto.ancestors(b).tolist()) == {b, onto.index["GO:0000001"], onto.index["GO:0008150"]}


GAF = "\t".join(["SGD", "S1", "ABC1", "", "GO:0000001", "PMID:1", "IDA", "", "P", "x", "YAL001C|abc1", "gene",
                 "taxon:4932", "2010", "SGD"]) + "\n" + \
      "\t".join(["SGD", "S2", "DEF1", "NOT", "GO:0000002", "PMID:2", "IMP", "", "P", "y", "YAL002W", "gene",
                 "taxon:4932", "2010", "SGD"]) + "\n" + \
      "\t".join(["SGD", "S3", "GHI1", "", "GO:0000002", "PMID:3", "IEA", "", "P", "z", "YAL003W|abc1", "gene",
                 "taxon:4932", "2010", "SGD"]) + "\n"


def test_gaf_and_mapper(tmp_path):
    p = tmp_path / "g.gaf.bz2"
    with bz2.open(p, "wt") as fh:
        fh.write("!gaf-version: 1.0\n" + GAF)
    recs = read_gaf(p)
    assert [r.object_id for r in recs] == ["S1"]  # NOT and IEA dropped
    allr = read_gaf(p, exclude_evidence=(), exclude_not=False)
    m = IdMapper.from_records(allr)
    assert m.map("yal001c") == "SGD:S1"
    assert m.map("DEF1") == "SGD:S2"
    assert m.map("abc1") == "SGD:S1"  # symbol level beats the ambiguous synonym
    assert m.map("nope") is None
    assert m.map_all(["N/A", "YAL003W"]) == ["SGD:S3"]


BIOGRID_V1 = """header junk
INTERACTOR_A\tINTERACTOR_B\tOFFICIAL_SYMBOL_A\tOFFICIAL_SYMBOL_B\tALIASES_FOR_A\tALIASES_FOR_B\tEXPERIMENTAL_SYSTEM\tSOURCE\tPUBMED_ID\tORGANISM_A_ID\tORGANISM_B_ID
YAL001C\tYAL002W\tABC1\tDEF1\tN/A\tN/A\tTwo-hybrid\tX\t1\t4932\t4932
YAL001C\tYAL003W\tABC1\tGHI1\tN/A\tN/A\tSynthetic Lethality\tX\t2\t4932\t4932
"""

BIOGRID_TAB3 = "#BioGRID Interaction ID\tEntrez Gene Interactor A\tEntrez Gene Interactor B\tBioGRID ID Interactor A\t" \
    "BioGRID ID Interactor B\tSystematic Name Interactor A\tSystematic Name Interactor B\tOfficial Symbol Interactor A\t" \
    "Official Symbol Interactor B\tSynonyms Interactor A\tSynonyms Interactor B\tExperimental System\t" \
    "Experimental System Type\tAuthor\tPublication Source\tOrganism ID Interactor A\tOrganism ID Interactor B\n" \
    "1\t11\t12\t1\t2\tYAL001C\tYAL002W\tABC1\tDEF1\t-\t-\tTwo-hybrid\tphysical\tX\tPUBMED:5\t559292\t559292\n" \
    "2\t11\t13\t1\t3\tYAL001C\tYAL003W\tABC1\tGHI1\t-\t-\tSynthetic Lethality\tgenetic\tX\tPUBMED:6\t559292\t559292\n"


def test_biogrid_formats(tmp_path):
    p1 = tmp_path / "b1.txt"
    p1.write_text(BIOGRID_V1)
    p3 = tmp_path / "b3.txt"
    p3.write_text(BIOGRID_TAB3)
    for p in (p1, p3):
        phys = read_biogrid(p)
        assert len(phys) == 1 and phys[0].names_a[0] == "YAL001C" and phys[0].names_b[0] == "YAL002W"
        assert len(read_biogrid(p, physical_only=False)) == 2
    assert len(read_biogrid(p3, pubmeds=["5"])) == 1
