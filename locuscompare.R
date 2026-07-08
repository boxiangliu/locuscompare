# NOTE: The plotting and LD/position helpers that used to live here have been
# replaced by the `locuscomparer` package (v1.1.0), which the app now depends on
# (see global.R). The Shiny server calls, from the package:
#   - get_position()      SNP coordinates (colotool tkg_p3v5a_hg19/hg38)
#   - retrieve_LD()       pairwise r^2 (colotool tkg_p3v5a_ld_chr*_*)
#   - get_lead_snp()      lead SNP selection
#   - assign_color()      LD-based point coloring
#   - make_scatterplot()  the LocusCompare scatter
#   - make_locuszoom()    the LocusZoom panels
# This file is intentionally left as a pointer and is no longer sourced.
