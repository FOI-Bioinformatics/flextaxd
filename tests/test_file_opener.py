import gzip

from flextaxd.modules import CreateCentrifugeDB, CreateGanonDB, CreateKrakenDatabase, functions


def test_the_opener_is_the_shared_one():
    """The gz sniffing opener existed four times, three of them unable to hand out bytes: the
    builders and the taxonomy readers now share one implementation (issue #117)"""
    for module in (CreateKrakenDatabase, CreateGanonDB, CreateCentrifugeDB):
        assert module.zopen is functions.zopen


def test_the_opener_returns_what_the_mode_asks_for(tmp_path):
    """A mode containing b stays binary through the gz sniffing, the other modes keep returning
    text, and gzipped files behave like unpacked ones (issue #117)"""
    text = "first line\nsecond line\n"
    plain = tmp_path / "plain.fna"
    with open(str(plain), "w") as handle:
        handle.write(text)
    packed = tmp_path / "packed.fna.gz"
    with gzip.open(str(packed), "wt") as handle:
        handle.write(text)
    for path in (str(plain), str(packed)):
        with functions.zopen(path, "rb") as handle:
            assert handle.read() == text.encode()
        with functions.zopen(path, "r") as handle:
            assert handle.read() == text
        with functions.zopen(path) as handle:
            assert handle.read() == text


def test_the_gz_suffix_is_removed_as_a_suffix():
    """.gz is a suffix and not a set of characters: rstrip(".gz") also removed the g of a name
    that merely ended in g, z or a point (issue #118)"""
    assert functions.strip_gz("GCF_000005845.2_ASM584v2_genomic.fna.gz") == "GCF_000005845.2_ASM584v2_genomic.fna"
    assert functions.strip_gz("batch_gz") == "batch_gz"
    assert functions.strip_gz("nt.fna") == "nt.fna"
