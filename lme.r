library(openxlsx)
library(stringr)
library(purrr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(patchwork)
library(forcats)
library(Seurat)
library(glmGamPoi)
library(transformGamPoi)
library(ggrepel)
library(data.table)

DT <- as.data.table(read.xlsx("data/Data_240411_initial_file.xlsx"))[
  , .(
    id = sprintf("%s_%s", Patient, Timepoint),
    patient = Patient,
    timepoint = Timepoint,
    t1.gd.num = as.numeric(str_trim(`T1.Gd#`)),
    t1.gd.vol = as.numeric(str_trim(T1Gd.Vol)),
    log.snfl = log10(as.numeric(str_trim(sNfL))),
    whole.brain.vol = as.numeric(str_trim(Whole.brain.vol)),
    grey.matter.vol = as.numeric(str_trim(Grey.matter.vol)),
    white.matter.vol = as.numeric(str_trim(White.matter.vol)),
    thalamic.vol = as.numeric(str_trim(Thal.vlo)),
    EDSS = as.numeric(str_trim(EDSS.Score))
  )
][order(patient, timepoint)][
  fread("data/ebv_load.tsv"),
  on = .(patient, timepoint),
  nomatch = NULL
][, `:=`(
  EBV.load = as.numeric(EBV.load),
  EBV.lytic = as.numeric(EBV.lytic),
  EBV.latent = as.numeric(EBV.latent)
)][
  order(timepoint),
  `:=`(
    t1.gd.num.next = shift(t1.gd.num, -1),
    t1.gd.vol.next = shift(t1.gd.vol, -1),
    log.snfl.next = shift(log.snfl, -1),
    timepoint = as.character(timepoint)
  ),
  by = .(patient)
][startsWith(patient, "OT") & timepoint %in% c("1", "3", "4", "6", "7")]



if (TRUE) {
  mem <- function(resp.var, data) {
    lmer(
      sprintf(
        # uncorrelated EBV load / latent / lytic
        "%s ~ EBV.load + (EBV.load | patient) + EBV.lytic + (EBV.lytic | patient) + EBV.latent + (EBV.latent | patient) + (1 | timepoint)",
        # correlated EBV load / latent / lytic
        # "%s ~ EBV.load * EBV.lytic * EBV.latent + ((EBV.load + EBV.lytic + EBV.latent) | patient) + (1 | timepoint)",
        resp.var
      ),
      data
    )
  }
  mod.cur <- list(
    t1.gd.num = mem("t1.gd.num", DT),
    t1.gd.vol = mem("t1.gd.vol", DT),
    snfl = mem("log.snfl", DT)
  )
  mod.fut <- list(
    t1.gd.num.next = mem("t1.gd.num.next", DT),
    t1.gd.vol.next = mem("t1.gd.vol.next", DT),
    snfl.next = mem("log.snfl.next", DT)
  )
  fixed <- rbindlist(imap(append(mod.cur, mod.fut), function(x, n) {
    as.data.table(summary(x)[["coefficients"]], keep.rownames = TRUE)[
      ,
      .(
        response.var = n,
        fixed.eff = rn,
        estimate = Estimate,
        p.val = `Pr(>|t|)`,
        signif = `Pr(>|t|)` < 0.05
      )
    ]
  }))
  random <- rbindlist(imap(append(mod.cur, mod.fut), function(x, resp.var) {
    imap(summary(x)[["varcor"]], function(x, random.eff) {
      x <- attr(x, "stddev")
      data.table(response.var = resp.var, random.eff = random.eff, intercept.var = x[[1]], slope.var = c(x, NA)[[2]], slope.for = c(names(x), NA)[[2]])
    })
  }) |> list_flatten())

  mms <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }

  lme.pl <- melt(
    DT[, .(
      id,
      .SD[, lapply(.SD, mms), .SDcols = patterns("EBV")]
    )],
    id.vars = "id",
    variable.name = "transcript.type",
    value.name = "level"
  )[melt(
    DT[, .(id, .SD[, lapply(.SD, mms), .SDcols = patterns("t1.gd|snfl")])],
    id.vars = "id",
    variable.name = "response.var",
    value.name = "response"
  ), on = .(id), allow.cartesian = TRUE] |>
    ggplot() +
    aes(x = response, y = level, colour = transcript.type) +
    geom_point() +
    facet_wrap(vars(response.var)) +
    labs(y = "scaled transcript level", x = "scaled response", colour = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

if (FALSE) {
  so <- readRDS("data/thebault_proc.seurat.rds")
  so <- so[, endsWith(so$cell.type, "B")]
  bulk.mtx <- as.matrix(AggregateExpression(so, group.by = c("patient", "timepoint"), assays = "RNA")[["RNA"]])
  fwrite(as.data.table(bulk.mtx, keep.rownames = "SYMBOL"), "data/b_cell_bulk.tsv", sep = "\t")
}

if (FALSE) {
  bulk.mtx <- as.matrix(fread("data/b_cell_bulk.tsv"), rownames = "SYMBOL")

  # drop samples where EBV transcript values were not recorded
  bulk.mtx <- bulk.mtx[, DT[!is.na(EBV.load), id]]

  if (FALSE) {
    # drop genes not present in at least 90% of samples
    min.expr <- ncol(bulk.mtx) * 0
    bulk.mtx <- bulk.mtx[apply(bulk.mtx, 1, function(x) {
      sum(x > 0) >= min.expr
    }), ]
  }

  fit <- glm_gp(
    data = bulk.mtx,
    design = as.formula("~ EBV.load + EBV.lytic + EBV.latent"),
    col_data = DT[!is.na(EBV.load), .(id, EBV.load, EBV.latent, EBV.lytic)] |>
      data.frame(row.names = "id") |>
      scale() |>
      as.data.frame(),
    size_factors = "ratio"
  )

  out.res <- map(c("EBV.load", "EBV.lytic", "EBV.latent"), function(x) {
    as.data.table(test_de(fit, x))[, `:=`(readout = x)][]
  }) |> rbindlist()
  fwrite(out.res, "glm_gp_out.tsv", sep = "\t")
}

if (FALSE) {
  degs.pl <- out.res[adj_pval < 0.05 & abs(lfc) < 10] |>
    ggplot() +
    aes(
      y = fct_reorder(name, abs(lfc)),
      x = lfc,
      colour = -log10(adj_pval),
      shape = fct_relevel(readout, "EBV.load", "EBV.lytic", "EBV.latent")
    ) +
    geom_point(size = 5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_colour_viridis_c(option = "rocket") +
    labs(shape = NULL, y = NULL)

  vlc.pl <- ggplot(out.res) +
    aes(
      x = lfc,
      y = -log10(adj_pval),
      label = fifelse(adj_pval < 0.05 & abs(lfc) > 1.5, name, NA),
      colour = fct_relevel(
        fcase(
          adj_pval >= 0.05, "ns",
          lfc > 0, "up",
          lfc < 0, "down"
        ),
        "up", "down", "ns"
      )
    ) +
    geom_point() +
    geom_label_repel(show.legend = FALSE) +
    labs(colour = NULL) +
    facet_wrap(vars(readout))

  expr.pl <- melt(
    as.data.table(colData(fit$data), keep.rownames = "id"),
    id.vars = "id",
    variable.name = "response.var",
    value.name = "scaled.resp"
  )[
    melt(
      as.data.table(
        t(transformGamPoi(fit, "shifted_log")),
        keep.rownames = "id"
      ),
      id.vars = "id",
      variable.name = "gene",
      value.name = "shifted_log.expr"
    ),
    on = .(id),
    allow.cartesian = TRUE
  ][
    out.res[adj_pval < 0.05 & abs(lfc) < 10],
    on = .(gene = name, response.var = readout)
  ] |>
    ggplot() +
    aes(x = shifted_log.expr, y = scaled.resp, colour = response.var) +
    geom_point() +
    facet_wrap(vars(gene, lfc)) +
    labs(
      x = "gene expression (shifted log)",
      y = "scaled readout value",
      colour = NULL
    )
}
