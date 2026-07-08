# Configuration for the LocusCompare single-locus Shiny app.
# Environment-driven: NO secrets committed to the repo. Provide values via the
# environment (see .env.example) — e.g. the container passes them through .env.
#
# The database is LocusCompare2's shared reference server (`colotool`): it holds the
# 1000 Genomes variant/LD reference tables (tkg_p3v5a_*) and hg38 gencode that the
# single-locus plots need. Read-only user. The password is the one bundled in the
# public `locuscomparer` package (data/config.rda); we still keep it out of git.
# Boxiang Liu

# Application directory (where ui.R/server.R/init.js/data live). shiny-server sets the
# working directory to the app dir; APP_HOME lets the container override it explicitly.
home_dir = Sys.getenv('APP_HOME', getwd())
data_dir = file.path(home_dir, 'data')

# Database connection (LocusCompare2 / colotool reference tables):
db_host = Sys.getenv('DB_HOST', '54.254.162.217')
db_port = as.integer(Sys.getenv('DB_PORT', '31987'))
db_name = Sys.getenv('DB_NAME', 'colotool')
db_user = Sys.getenv('DB_USER', 'locuscompare')
db_pass = Sys.getenv('DB_PASSWORD', '')
