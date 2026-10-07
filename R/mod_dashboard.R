mod_dashboard_ui <- function(id, data) {
  ns <- NS(id)
  patients <- data$patients
  select_filter <- function(name, label, values, missing_label = "Not recorded") {
    choices <- revive_filter_choices(values, missing_label)
    selectizeInput(
      ns(name), label, choices = choices[choices != "__all__"],
      selected = character(0), multiple = TRUE,
      options = list(searchField = "label", placeholder = "All", plugins = list("remove_button"))
    )
  }
  layout <- bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      id = ns("sidebar"), width = 280, open = "desktop", bg = "#F4FCFB",
      tags$div(class = "filter-heading", tags$h2("Filters"), icon("sliders")),
      tags$p(class = "filter-intro", "Define the patients in your view."),
      select_filter("study_group", "Study group", patients$study_group),
      select_filter("infection", "Infection", patients$infection, "No infection reported"),
      select_filter("vaccination", "Vaccination", patients$vaccination, "No vaccination reported"),
      select_filter("sex", "Sex", patients$sex),
      select_filter("sotr", "SOTR", patients$sotr, "No transplant organs recorded"),
      tags$p(class = "filter-hint", "Organs from the most recent transplant"),
      select_filter("immunosuppression", "Immunosuppression", patients$immunosuppression),
      sliderInput(
        ns("age"), tags$span(class = "age-filter-label", "Age", textOutput(ns("age_label"), inline = TRUE)),
        min = data$age_limits[1], max = data$age_limits[2], value = data$age_limits,
        step = 1, sep = "", width = "100%"
      ),
      tags$p(class = "filter-hint age-hint", "Age at consent, in years. All includes ages not recorded."),
      actionButton(ns("reset"), "Reset filters", icon = icon("arrow-rotate-left"), class = "reset-filters"),
      tags$div(class = "sidebar-note", icon("circle-info"),
        tags$span("Select one or more categories per filter. Patients can match any selected category within a filter. Filters apply together; an empty filter includes all patients."))
    ),
    fillable = FALSE, fill = FALSE, border = FALSE, padding = 24, gap = 24,
    tags$main(class = "dashboard-main",
      tags$div(class = "dashboard-intro",
        tags$div(tags$p(class = "eyebrow", "PARTICIPANT EXPLORER"),
          tags$h2("Cohort overview"),
          tags$p("Explore patient characteristics and events over time.")),
        tags$span(class = "study-badge", tags$span(class = "status-dot"), "Longitudinal study")
      ),
      tags$div(class = "metric-grid",
        bslib::card(class = "metric-card metric-patients",
          tags$div(class = "metric-heading", "Patients", icon("users")),
          textOutput(ns("patient_count"), container = function(...) tags$div(class = "metric-value", ...)),
          tags$p(class = "metric-caption", paste("of", nrow(patients), "patients in the full cohort"))
        ),
        bslib::card(class = "metric-card metric-secondary",
          tags$div(class = "metric-heading", "Lorem ipsum", icon("chart-simple")),
          tags$div(class = "metric-placeholder", "Lorem ipsum"),
          tags$p(class = "metric-caption", "Lorem ipsum dolor sit amet, consectetur adipiscing elit.")
        ),
        bslib::card(class = "metric-card metric-tertiary",
          tags$div(class = "metric-heading", "Lorem ipsum", icon("layer-group")),
          tags$div(class = "metric-placeholder", "Lorem ipsum"),
          tags$p(class = "metric-caption", "Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.")
        )
      ),
      bslib::card(class = "patient-card",
        bslib::card_header(
          tags$div(tags$h3("Patient details"), tags$p("One row per patient in the selected cohort.")),
          actionButton(ns("export_options"), "Export CSV", icon = icon("download"), class = "export-button")
        ),
        bslib::card_body(DT::DTOutput(ns("patients_table")))
      ),
      bslib::card(class = "timeline-card",
        bslib::card_header(
          tags$div(tags$h3("Infection & vaccination timeline"),
            tags$p("All recorded events for the selected patients.")),
          tags$div(class = "event-legend", `aria-label` = "Event legend",
            tags$span(tags$i(class = "legend-infection", `aria-hidden` = "true"), "Infection"),
            tags$span(tags$i(class = "legend-vaccination", `aria-hidden` = "true"), "Vaccination"))
        ),
        bslib::card_body(
          tags$div(class = "timeline-summary", textOutput(ns("timeline_summary"))),
          tags$div(class = "timeline-scroll", uiOutput(ns("swimmer_container"))),
          tags$p(class = "timeline-note", "Lines connect first and last recorded events; follow-up duration is not available. Reports on the same date may overlap.")
        )
      ),
      tags$footer(class = "dashboard-footer", "REVIVE", tags$span("Patient and event explorer"))
    )
  )
  tags$div(class = "revive-dashboard", layout)
}

mod_dashboard_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    filter_names <- c("study_group", "infection", "vaccination", "sex", "sotr", "immunosuppression")
    export_columns <- revive_export_columns()
    filtered_patients <- reactive({
      filters <- stats::setNames(lapply(filter_names, function(name) input[[name]]), filter_names)
      filter_revive_patients(data$patients, filters, input$age, data$age_limits)
    })
    patient_table <- reactive(revive_patient_table(filtered_patients()))
    filtered_events <- reactive({
      data$events[data$events$record_id %in% filtered_patients()$record_id, , drop = FALSE]
    })
    dated_events <- reactive({
      events <- filtered_events()
      events[!is.na(events$event_date), , drop = FALSE]
    })

    output$patient_count <- renderText(format(nrow(filtered_patients()), big.mark = ","))
    output$age_label <- renderText({
      if (is.null(input$age) || identical(as.numeric(input$age), as.numeric(data$age_limits))) "All"
      else paste(input$age, collapse = "–")
    })
    observeEvent(input$reset, {
      for (name in filter_names) updateSelectizeInput(session, name, selected = character(0))
      updateSliderInput(session, "age", value = data$age_limits)
    })

    output$patients_table <- DT::renderDT({
      DT::datatable(
        patient_table(), rownames = FALSE, selection = "none", escape = TRUE,
        style = "bootstrap5", class = "table table-hover", width = "100%",
        options = list(
          pageLength = 5, lengthChange = FALSE, searching = FALSE, scrollX = TRUE,
          dom = "t<'table-bottom'ip>", autoWidth = FALSE, order = list(list(0, "asc")),
          columnDefs = list(list(targets = 0, type = "string")),
          language = list(
            emptyTable = "No patients match these filters. Try adjusting or resetting the filters.",
            info = "Showing _START_–_END_ of _TOTAL_ patients", infoEmpty = "0 patients",
            paginate = list(previous = "Previous", `next` = "Next")
          ),
          createdRow = DT::JS("function(row, data) { $('td', row).each(function(i) { this.title = data[i] == null ? '' : data[i]; }); }")
        )
      )
    }, server = TRUE)

    observeEvent(input$export_options, {
      selected_columns <- intersect(input$export_columns, export_columns)
      showModal(modalDialog(
        title = "Export patient data",
        tags$p(paste(format(nrow(filtered_patients()), big.mark = ","),
                     "patients match the current filters.")),
        tags$p("All displayed columns are included. Choose any additional variables to append to the CSV."),
        selectizeInput(
          session$ns("export_columns"), "Additional variables",
          choices = NULL, multiple = TRUE,
          options = list(placeholder = "Search variables", maxOptions = 100,
                         plugins = list("remove_button"))
        ),
        tags$p("Leave the selection empty to download only the displayed columns."),
        footer = tagList(
          modalButton("Cancel"),
          downloadButton(session$ns("export_csv"), "Download CSV", class = "export-button")
        ),
        easyClose = TRUE
      ))
      updateSelectizeInput(
        session, "export_columns", choices = export_columns,
        selected = selected_columns, server = TRUE
      )
    })

    output$export_csv <- downloadHandler(
      filename = function() paste0("REVIVE_patients_", Sys.Date(), ".csv"),
      content = function(file) {
        export <- revive_patient_export(
          filtered_patients(), data$patient_characteristics, input$export_columns
        )
        utils::write.csv(export, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
      },
      contentType = "text/csv; charset=UTF-8"
    )
    output$timeline_summary <- renderText({
      events <- dated_events()
      missing_dates <- nrow(filtered_events()) - nrow(events)
      paste0(
        length(unique(events$record_id)), " of ", nrow(filtered_patients()),
        " selected patients have dated events · ", nrow(events), " event reports",
        if (missing_dates) paste0(" · ", missing_dates, " undated reports omitted") else ""
      )
    })
    output$swimmer_container <- renderUI({
      if (!nrow(dated_events())) {
        return(tags$div(class = "timeline-empty", role = "status",
          icon("calendar", `aria-hidden` = "true"),
          tags$strong("No dated events for the selected patients."),
          tags$p("Adjust the filters to explore another group of patients.")
        ))
      }
      height <- max(340, length(unique(dated_events()$record_id)) * 25 + 100)
      plotly::plotlyOutput(session$ns("swimmer_plot"), height = paste0(height, "px"))
    })
    output$swimmer_plot <- plotly::renderPlotly({
      req(nrow(dated_events()) > 0L)
      revive_swimmer_plot(dated_events())
    })
    list(patients = filtered_patients, table = patient_table, events = filtered_events)
  })
}
