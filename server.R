# Server for the LocusCompare single-locus app.
# Plotting, lead-SNP selection, and LD coloring come from the `locuscomparer`
# package (pure functions, no DB access). ALL database access — region variants,
# LD, gene coordinates, SNP info — goes through the app's connection pool to
# `colotool` (config.R / env), so there is a single source of DB configuration.
# Boxiang Liu

options(shiny.maxRequestSize = 5 * 1024^2)

# ---- helpers ---------------------------------------------------------------

# Read an uploaded association file; must have columns rsid, pval.
read_upload <- function(datapath) {
    d <- data.table::fread(datapath, header = TRUE)
    if (!all(c('rsid', 'pval') %in% colnames(d))) return(NULL)
    data.table::data.table(rsid = as.character(d$rsid), pval = as.numeric(d$pval))
}

# Empty LD frame (so assign_color() degrades to all-blue when a SNP has no partners).
empty_ld <- function() data.frame(SNP_A = character(), SNP_B = character(), R2 = numeric())

shinyServer(function(input, output, session) {

    message('SESSION STARTED')

    # Plots tab is revealed only after a plot is requested.
    hideTab(inputId = 'navbarPage', target = 'Plots')

    tmp_dir <- file.path(tempdir(), session$token)
    dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

    observeEvent(input$go_to_website, {
        shinyjs::hide(id = "select-genome", anim = TRUE, animType = "fade")
        shinyjs::show(id = "app-content")
    })

    # ---- assembly / reference tables --------------------------------------
    build <- reactive({ if (input$genome_assembly == 'GRCh37/hg19') 'hg19' else 'hg38' })
    genome <- reactive({ paste0('tkg_p3v5a_', build()) })
    transcriptome <- reactive({
        if (build() == 'hg19') 'gencode_v19_gtex_v6p' else 'gencode_v19_gtex_v6p_hg38'
    })

    # ---- input validation --------------------------------------------------
    study1_ok <- reactive({ isTruthy(input$file1) & isTruthy(input$file1_trait) })
    study2_ok <- reactive({ isTruthy(input$file2) & isTruthy(input$file2_trait) })

    output$check_study1 <- renderText({
        if (study1_ok()) 'Study 1 has a valid input.' else 'Please upload a file and enter a label for Study 1.'
    })
    output$check_study2 <- renderText({
        if (study2_ok()) 'Study 2 has a valid input.' else 'Please upload a file and enter a label for Study 2.'
    })

    selected_snp_region <- reactive({ isTruthy(input$reference_snp) & isTruthy(input$snp_window) })
    selected_gene_region <- reactive({ isTruthy(input$reference_gene) & isTruthy(input$gene_window) })
    selected_coordinate <- reactive({ isTruthy(input$chr) & isTruthy(input$start) & isTruthy(input$end) })
    region_total <- reactive({ selected_snp_region() + selected_gene_region() + selected_coordinate() })

    output$check_region <- renderText({
        n <- region_total()
        if (n == 0) 'Please type in a SNP, a gene, or a coordinate.'
        else if (n == 1) 'Region has a valid input.'
        else 'Please select only one of SNP, gene, or coordinate.'
    })

    ready <- reactive({ study1_ok() & study2_ok() & region_total() == 1 })
    observe({ shinyjs::toggleState('interactive_to_locuscompare', ready()) })

    # ---- navigation --------------------------------------------------------
    observeEvent(input$interactive_to_locuscompare, {
        req(ready())
        showTab(inputId = 'navbarPage', target = 'Plots', select = TRUE)
    })
    observeEvent(input$back, {
        updateNavbarPage(session, 'navbarPage', selected = 'Single Locus')
        hideTab(inputId = 'navbarPage', target = 'Plots')
    })

    # ---- region ------------------------------------------------------------
    region <- eventReactive(input$interactive_to_locuscompare, {
        req(ready())

        if (selected_snp_region()) {
            chr_pos <- dbGetQuery(locuscompare_pool, sprintf(
                'select chr, pos from %s where rsid = "%s" limit 1;', genome(), input$reference_snp))
            shiny::validate(need(nrow(chr_pos) == 1, sprintf('SNP %s not found.', input$reference_snp)))
            w <- max(1, min(input$snp_window, 1000))
            list(chr = as.character(chr_pos$chr), start = chr_pos$pos - w * 500, end = chr_pos$pos + w * 500)

        } else if (selected_gene_region()) {
            gene <- tryCatch(
                dbGetQuery(locuscompare_pool, sprintf(
                    'select chr, start, end from %s where gene_name = "%s" limit 1;', transcriptome(), input$reference_gene)),
                error = function(e) NULL)
            shiny::validate(need(!is.null(gene),
                'Gene lookup is unavailable for hg19 on this server; use hg38 or select by SNP/coordinate.'))
            shiny::validate(need(nrow(gene) == 1, sprintf('Gene %s not found.', input$reference_gene)))
            w <- max(1, min(input$gene_window, 1000))
            mid <- (gene$start + gene$end) / 2
            list(chr = sub('^chr', '', as.character(gene$chr)), start = mid - w * 500, end = mid + w * 500)

        } else {
            start_ <- max(1, input$start)
            end_ <- max(start_, min(input$end, start_ + 1e6))
            list(chr = as.character(input$chr), start = start_, end = end_)
        }
    })

    chr <- reactive({ region()$chr })

    # Variants (rsid + position) in the selected window, from colotool.
    region_variants <- reactive({
        r <- region()
        v <- dbGetQuery(locuscompare_pool, sprintf(
            'select rsid, pos from %s where chr = "%s" and pos >= %d and pos <= %d;',
            genome(), r$chr, as.integer(r$start), as.integer(r$end)))
        data.table::as.data.table(v)
    })

    # ---- data --------------------------------------------------------------
    d1 <- reactive({
        req(input$file1)
        d <- read_upload(input$file1$datapath)
        shiny::validate(need(!is.null(d), 'Study 1 file must have columns: rsid, pval.'))
        d
    })
    d2 <- reactive({
        req(input$file2)
        d <- read_upload(input$file2$datapath)
        shiny::validate(need(!is.null(d), 'Study 2 file must have columns: rsid, pval.'))
        d
    })

    merged <- eventReactive(input$interactive_to_locuscompare, {
        req(ready())
        rv <- region_variants()
        shiny::validate(need(nrow(rv) > 0, 'No reference variants found in the selected region.'))

        m <- merge(as.data.frame(d1()), as.data.frame(d2()), by = 'rsid', suffixes = c('1', '2'))
        m <- merge(m, as.data.frame(rv), by = 'rsid')
        shiny::validate(need(nrow(m) > 0,
            'No overlapping SNPs between the two studies within the selected region.'))

        m$chr <- region()$chr
        m$logp1 <- -log10(m$pval1)
        m$logp2 <- -log10(m$pval2)
        m
    })

    # ---- lead SNP selection ------------------------------------------------
    snp <- reactiveVal('')

    observeEvent(merged(), {
        m <- merged()
        lead <- locuscomparer::get_lead_snp(m)
        snp(lead)
        zoom$xmin <- NULL; zoom$xmax <- NULL          # reset any prior zoom for the new region
        updateSelectizeInput(session, 'snp', choices = m$rsid, selected = lead, server = TRUE)
    })
    observeEvent(input$snp, { if (isTruthy(input$snp)) snp(input$snp) }, ignoreInit = TRUE)
    observeEvent(input$plot_click, {
        req(merged())
        picked <- nearPoints(plot_data(), input$plot_click, xvar = 'logp1', yvar = 'logp2',
                             maxpoints = 1, threshold = 10)$rsid
        if (length(picked) == 1) {
            snp(picked)
            updateSelectizeInput(session, 'snp', selected = picked)
        }
    })

    # ---- LD + aesthetics ---------------------------------------------------
    # LD via the app pool (single DB config), same tables locuscomparer would query.
    ld <- reactive({
        req(snp() != '', chr())
        cc <- chr(); ss <- snp(); pop <- input$population
        out <- tryCatch({
            r1 <- dbGetQuery(locuscompare_pool, sprintf(
                "select SNP_A, SNP_B, R2 from tkg_p3v5a_ld_chr%s_%s where SNP_A = '%s';", cc, pop, ss))
            r2 <- dbGetQuery(locuscompare_pool, sprintf(
                "select SNP_B as SNP_A, SNP_A as SNP_B, R2 from tkg_p3v5a_ld_chr%s_%s where SNP_B = '%s';", cc, pop, ss))
            rbind(r1, r2)
        }, error = function(e) NULL)
        if (is.null(out) || nrow(out) == 0) empty_ld() else out
    })

    color <- reactive({ locuscomparer::assign_color(merged()$rsid, snp(), ld()) })
    shape <- reactive({ s <- ifelse(merged()$rsid == snp(), 23, 21); names(s) <- merged()$rsid; s })
    size  <- reactive({ s <- ifelse(merged()$rsid == snp(), 3, 2);  names(s) <- merged()$rsid; s })

    # ---- zoom + plot data --------------------------------------------------
    zoom <- reactiveValues(xmin = NULL, xmax = NULL)
    observeEvent(input$plot_dblclick, {
        b <- input$plot_brush
        if (!is.null(b)) { zoom$xmin <- b$xmin; zoom$xmax <- b$xmax }
        else { zoom$xmin <- NULL; zoom$xmax <- NULL }
    })

    plot_data <- reactive({
        pd <- merged()
        if (!is.null(zoom$xmin)) pd <- pd[pd$pos >= zoom$xmin & pd$pos <= zoom$xmax, ]
        pd$label <- ifelse(pd$rsid == snp(), pd$rsid, '')
        pd
    })

    title1 <- reactive({ input$file1_trait })
    title2 <- reactive({ input$file2_trait })

    # ---- outputs -----------------------------------------------------------
    output$blank_plot <- renderPlot({ ggplot() + geom_blank() })

    locuscompare_plot <- reactive({
        pd <- plot_data()
        shiny::validate(need(nrow(pd) > 0, 'No SNPs to display in the current view.'))
        locuscomparer::make_scatterplot(pd, title1(), title2(), color(), shape(), size(), legend = FALSE)
    })
    locuszoom_plot <- function(y_col, title) {
        pd <- plot_data()
        pd$logp <- pd[[y_col]]
        locuscomparer::make_locuszoom(pd, title, chr(), color(), shape(), size())
    }

    output$locuscompare <- renderPlot({ locuscompare_plot() })
    output$locuszoom1 <- renderPlot({
        shiny::validate(need(nrow(plot_data()) > 0, 'No SNPs to display.'))
        locuszoom_plot('logp1', title1())
    })
    output$locuszoom2 <- renderPlot({
        shiny::validate(need(nrow(plot_data()) > 0, 'No SNPs to display.'))
        locuszoom_plot('logp2', title2())
    })

    output$snp_info <- renderText({
        req(snp() != '')
        res <- dbGetQuery(locuscompare_pool, sprintf(
            'select chr,pos,rsid,ref,alt,AF,AFR_AF,AMR_AF,EAS_AF,EUR_AF,SAS_AF from %s where rsid = "%s" limit 1;',
            genome(), snp()))
        if (nrow(res) == 0) return(sprintf('rs ID: %s (no reference record)', snp()))
        sprintf(paste0('Chromosome: %s\nPosition: %s\nrs ID: %s\nReference: %s\nAlternate: %s\n',
                       'Allele Frequency: %s\nAFR: %s\nAMR: %s\nEAS: %s\nEUR: %s\nSAS: %s'),
                res$chr, res$pos, res$rsid, res$ref, res$alt,
                res$AF, res$AFR_AF, res$AMR_AF, res$EAS_AF, res$EUR_AF, res$SAS_AF)
    })

    output$ld_snps <- renderDataTable({
        req(merged())
        ld_tab <- ld()
        ld_tab <- ld_tab[ld_tab$SNP_A == snp(), c('SNP_B', 'R2')]
        colnames(ld_tab) <- c('rsid', 'r2')
        ld_tab <- rbind(data.frame(rsid = snp(), r2 = 1), ld_tab)
        out <- merge(
            merged()[, c('rsid', 'chr', 'pos', 'pval1', 'pval2')],
            ld_tab, by = 'rsid')
        out <- out[order(-out$r2), ]
        colnames(out) <- c('rsID', 'Chromosome', 'Position', 'P (x-axis)', 'P (y-axis)', 'LD (r2)')
        datatable(out, selection = 'none', rownames = FALSE)
    })

    # ---- example files -----------------------------------------------------
    output$file1_example <- downloadHandler(
        filename = function() 'PHACTR1_Artery_Coronary.tsv',
        content = function(file) file.copy(file.path(home_dir, 'data/example/PHACTR1_Artery_Coronary.tsv'), file))
    output$file2_example <- downloadHandler(
        filename = function() 'PHACTR1_Coronary_Heart_Disease_Nikpay_2015.tsv',
        content = function(file) file.copy(file.path(home_dir, 'data/example/PHACTR1_Coronary_Heart_Disease_Nikpay_2015.tsv'), file))

    # ---- download results --------------------------------------------------
    output$single_download <- downloadHandler(
        filename = function() 'results.zip',
        content = function(file) {
            data_fn <- file.path(tmp_dir, 'data.tsv')
            ld_fn <- file.path(tmp_dir, 'ld.tsv')
            lc_fn <- file.path(tmp_dir, 'locuscompare.pdf')
            lz1_fn <- file.path(tmp_dir, 'locuszoom1.pdf')
            lz2_fn <- file.path(tmp_dir, 'locuszoom2.pdf')

            data.table::fwrite(merged(), data_fn, sep = '\t')
            data.table::fwrite(ld(), ld_fn, sep = '\t')

            ggsave(lc_fn, locuscompare_plot(), width = input$locuscompare_length, height = input$locuscompare_length)
            ggsave(lz1_fn, locuszoom_plot('logp1', title1()), width = input$locuszoom_width, height = input$locuszoom_height)
            ggsave(lz2_fn, locuszoom_plot('logp2', title2()), width = input$locuszoom_width, height = input$locuszoom_height)

            utils::zip(file, c(data_fn, ld_fn, lc_fn, lz1_fn, lz2_fn), flags = '-j')
        })
})
