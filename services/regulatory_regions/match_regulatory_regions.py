def match_regulatory_regions(db, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a 6-long list where each regulatory region is set to 1 if the sample within it and 0 if not.
    The classifications order:  ['CTCF_binding_site', 'promoter', 'promoter_flanking_region', 'open_chromatin_region',
                                'TF_binding_site', 'enhancer']
    """
    res = db.get(f"{chromosome}-{start_pos}", 0) | db.get(f"{chromosome}-{end_pos}", 0)

    return [(res >> i) & 1 for i in range(5, -1, -1)]


