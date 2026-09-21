ui <- bslib::page_fluid(
  title = "REVIVE",
  lang = "en",
  theme = bslib::bs_theme(
    version = 5, bg = "#F4FCFB", fg = "#243C4C", primary = "#5289AD",
    secondary = "#698696", base_font = "Arial, sans-serif",
    heading_font = "Arial, sans-serif", border_radius = "0.8rem"
  ),
  tags$head(tags$link(rel = "stylesheet", href = "custom.css")),
  tags$header(
    class = "revive-header",
    tags$div(class = "revive-brand",
      tags$div(class = "revive-mark", `aria-hidden` = "true", icon("wave-square")),
      tags$div(
        tags$h1("REVIVE"),
        tags$p(tags$strong("Respiratory Viral Infection Vaccine Effectiveness in Transplant Recipients"))
      )
    ),
    tags$span(class = "revive-header-badge", "COHORT DASHBOARD")
  ),
  mod_dashboard_ui("dashboard", revive_data)
)
