
library(shiny)
library(dplyr)
library(readr)
library(igraph)
library(plotly)
library(tibble)


ui <- fluidPage(
  titlePanel('Forschungslandschaft Schweiz'),
  p('Forschende nach Forschungsfeld und institutioneller Zugehörigkeit'),
  selectInput(
    "gruppe",
    "Bitte wählen Sie Ihre Gruppe:",
    choices = c(
      "Gruppe 1" = "0",
      "Gruppe 2" = "1831865388",
      "Gruppe 3" = "508919860"
    )
  ),
  plotlyOutput('network', height = '850px'),
  h4('Aktuelle Daten'),
  tableOutput('data')
#  ),
)

# Repel label rectangles in data coordinates, including field headings.
# Approximate label widths depend on character counts; no extra packages needed.
repel_labels <- function(labels, xlim, ylim, iterations = 350L) {
  if (!nrow(labels)) return(labels)
  span_x <- diff(xlim)
  span_y <- diff(ylim)
  # Coordinates are normalized to make collision checks independent of plot aspect.
  ax <- (labels$x - xlim[1]) / span_x
  ay <- (labels$y - ylim[1]) / span_y
  lx <- ax
  ly <- ay + ifelse(labels$kind == 'Field', 0.025, 0.014)
  widths <- pmin(0.27, pmax(0.022, nchar(labels$label, type = 'width') *
                             ifelse(labels$kind == 'Field', 0.0095, 0.007)))
  heights <- ifelse(labels$kind == 'Field', 0.031, 0.024)
  n <- nrow(labels)
  for (iter in seq_len(iterations)) {
    dx <- numeric(n)
    dy <- numeric(n)
    for (i in seq_len(n)) {
      if (i == n) next
      for (j in (i + 1L):n) {
        ox <- (widths[i] + widths[j]) / 2 + 0.007 - abs(lx[i] - lx[j])
        oy <- (heights[i] + heights[j]) / 2 + 0.005 - abs(ly[i] - ly[j])
        if (ox > 0 && oy > 0) {
          if (ox < oy * 1.7) {
            dir <- if (lx[i] <= lx[j]) -1 else 1
            push <- min(ox, 0.018) * 0.55
            dx[i] <- dx[i] + dir * push
            dx[j] <- dx[j] - dir * push
          } else {
            dir <- if (ly[i] <= ly[j]) -1 else 1
            push <- min(oy, 0.014) * 0.55
            dy[i] <- dy[i] + dir * push
            dy[j] <- dy[j] - dir * push
          }
        }
      }
    }
    # Keep labels close to their own nodes, but allow them to move.
    lx <- lx + dx + 0.015 * (ax - lx)
    ly <- ly + dy + 0.015 * (ay + ifelse(labels$kind == 'Field', 0.025, 0.014) - ly)
    lx <- pmax(widths / 2 + 0.005, pmin(1 - widths / 2 - 0.005, lx))
    ly <- pmax(heights / 2 + 0.005, pmin(1 - heights / 2 - 0.005, ly))
  }
  labels$lx <- xlim[1] + lx * span_x
  labels$ly <- ylim[1] + ly * span_y
  labels
}

server <- function(input, output, session) {
  
  data <- reactive({
    invalidateLater(60000, session)
    
    raw <- tryCatch({
      
      sheet_url <- paste0(
        "https://docs.google.com/spreadsheets/d/",
        "1y3qQ8eUz6yw8-pYHxwp87Dqy6Y3_cIdc6ZmbXjWCNcA",
        "/export?format=csv&gid=",
        input$gruppe
      )
      
      read_csv(sheet_url, show_col_types = FALSE)
      
    }, error = function(e) NULL)
    
    validate(need(
      !is.null(raw),
      "Google Sheet konnte nicht geladen werden. Bitte Freigabe und CSV-Link prüfen."
    ))
    
    names(raw) <- tolower(trimws(names(raw)))
    
    required <- c(
      "institution", "person",
      "schwerpunkt", "forschungsfeld"
    )
    
    validate(need(
      all(required %in% names(raw)),
      paste(
        "Fehlende Spalten:",
        paste(setdiff(required, names(raw)), collapse = ", ")
      )
    ))

    raw |>
      mutate(across(all_of(required), ~ trimws(as.character(.x)))) |>
      mutate(
        person = na_if(person, ""),
        person = if_else(is.na(person), "Unbekannt", person),
        institution = if_else(
          is.na(institution) | institution == "",
          "Unbekannt", institution
        ),
        schwerpunkt = if_else(
          is.na(schwerpunkt) | schwerpunkt == "",
          "Nicht angegeben", schwerpunkt
        ),
        forschungsfeld = if_else(
          is.na(forschungsfeld) | forschungsfeld == "",
          "Nicht angegeben", forschungsfeld
        )
      )
  })

  network_data <- reactive({
    df <- data()
    validate(need(nrow(df) > 0, 'Keine gültigen Personen- und Forschungsfelddaten vorhanden.'))

    # Each person must belong to one field. Conflicts are resolved deterministically
    # by selecting the most frequent recorded field (alphabetical tie-break).
    persons <- df |>
      count(person, forschungsfeld, name = 'frequency') |>
      arrange(person, desc(frequency), forschungsfeld) |>
      group_by(person) |>
      slice_head(n = 1) |>
      ungroup() |>
      transmute(person, field = forschungsfeld, name = paste0('p:', person))

    affiliations <- df |>
      distinct(person, institution) |>
      inner_join(persons |> select(person), by = 'person')

    institutions <- affiliations |>
      distinct(institution) |>
      mutate(name = paste0('i:', institution))

    vertices <- bind_rows(
      persons |> transmute(name, label = person, field, type = FALSE),
      institutions |> transmute(name, label = institution,
                                field = NA_character_, type = TRUE)
    )
    edge_table <- affiliations |>
      transmute(from = paste0('p:', person), to = paste0('i:', institution)) |>
      distinct()
    g <- graph_from_data_frame(edge_table, directed = FALSE, vertices = vertices)
    stopifnot(is_bipartite(g))

    fields <- sort(unique(persons$field))
    n_fields <- length(fields)
    field_counts <- persons |> count(field, name = 'n') |>
      mutate(radius = pmax(3.5, 1.5 * sqrt(n) + 1.5))

    # Distribute cluster centers around a ring; handle a single field separately.
    if (n_fields == 1L) {
      centers <- field_counts |> mutate(cx = 0, cy = 0)
    } else {
      angles <- 2 * pi * (seq_len(n_fields) - 1) / n_fields
      ring_radius <- max(12, max(field_counts$radius) * n_fields / pi * 1.5)
      centers <- field_counts |>
        mutate(cx = ring_radius * cos(angles), cy = ring_radius * sin(angles))
    }

    # Stable sunflower layout for people in each field.
    person_xy <- persons |>
      left_join(centers, by = 'field') |>
      group_by(field) |>
      arrange(person, .by_group = TRUE) |>
      mutate(
        k = row_number(), total = n(),
        angle = k * pi * (3 - sqrt(5)),
        r = if_else(total == 1L, 0,
                    0.76 * radius * sqrt((k - 1) / pmax(total - 1, 1))),
        x = cx + r * cos(angle),
        y = cy + r * sin(angle)
      ) |>
      ungroup()

    # Institutions are located at the centroids of their affiliated people.
    institution_xy <- affiliations |>
      left_join(person_xy |> select(person, x, y), by = 'person') |>
      group_by(institution) |>
      summarise(x = mean(x), y = mean(y), n_people = n(), .groups = 'drop') |>
      mutate(name = paste0('i:', institution))

    positions <- bind_rows(
      person_xy |> select(name, x, y),
      institution_xy |> select(name, x, y)
    )
    nodes <- vertices |> left_join(positions, by = 'name')

    person_details <- df |>
      group_by(person) |>
      summarise(
        schwerpunkt = paste(sort(unique(schwerpunkt)), collapse = ', '),
        institution = paste(sort(unique(institution)), collapse = ', '),
        .groups = 'drop'
      )
    nodes <- nodes |>
      left_join(person_details, by = c('label' = 'person')) |>
      mutate(
        tooltip = if_else(
          !type,
          paste0('<b>', label, '</b><br>Forschungsfeld: ', field,
                 '<br>Institution: ', institution,
                 '<br>Schwerpunkt: ', schwerpunkt),
          paste0('<b>Institution: ', label, '</b>')
        )
      )

    edge_xy <- edge_table |>
      left_join(nodes |> select(name, x, y), by = c('from' = 'name')) |>
      rename(x0 = x, y0 = y) |>
      left_join(nodes |> select(name, x, y), by = c('to' = 'name')) |>
      rename(x1 = x, y1 = y)

    # Irregular convex hulls, padded around persons (rather than circles).
    # Works for clusters with one or two people as well.
    hulls <- bind_rows(lapply(split(person_xy, person_xy$field), function(sub) {
      ang <- seq(0, 2 * pi, length.out = 25)[-25]
      pad <- 1.6
      expanded <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
        data.frame(x = sub$x[i] + pad * cos(ang),
                   y = sub$y[i] + pad * sin(ang))
      }))
      idx <- chull(expanded$x, expanded$y)
      idx <- c(idx, idx[1])
      tibble(field = sub$field[1], x = expanded$x[idx], y = expanded$y[idx])
    }))

    palette <- grDevices::hcl.colors(n_fields, palette = 'Dark 3')
    colors <- setNames(palette, fields)
    list(nodes = nodes, edges = edge_xy, centers = centers,
         hulls = hulls, colors = colors, graph = g)
  })


  output$network <- renderPlotly({
    net <- network_data()
    nodes <- net$nodes
    colors <- net$colors
    p <- plot_ly()

    # Filled irregular community polygons, behind all edges and nodes.
    for (f in names(colors)) {
      h <- net$hulls |> filter(field == f)
      p <- p |> add_trace(
        x = h$x, y = h$y, type = 'scatter', mode = 'lines',
        fill = 'toself',
        fillcolor = grDevices::adjustcolor(colors[[f]], alpha.f = 0.12),
        line = list(color = grDevices::adjustcolor(colors[[f]], alpha.f = 0.55), width = 1.4),
        hoverinfo = 'skip', showlegend = FALSE, inherit = FALSE
      )
    }

    # Affiliation ties are permanently visible (no opacity slider).
    e <- net$edges
    if (nrow(e) > 0) {
      p <- p |> add_trace(
        x = as.vector(rbind(e$x0, e$x1, NA_real_)),
        y = as.vector(rbind(e$y0, e$y1, NA_real_)),
        type = 'scatter', mode = 'lines',
        line = list(color = 'rgba(75,85,95,0.55)', width = 1.3),
        hoverinfo = 'skip', showlegend = FALSE, inherit = FALSE
      )
    }

    # Nodes are plotted without attached text: labels are repelled separately.
    for (f in names(colors)) {
      sub <- nodes |> filter(!type, field == f)
      p <- p |> add_trace(
        data = sub, x = ~x, y = ~y, type = 'scatter', mode = 'markers',
        hovertext = ~tooltip, hoverinfo = 'text',
        marker = list(color = colors[[f]], size = 9,
                      line = list(color = 'white', width = 1)),
        name = f, showlegend = TRUE, inherit = FALSE
      )
    }
    inst <- nodes |> filter(type)
    if (nrow(inst) > 0) {
      p <- p |> add_trace(
        data = inst, x = ~x, y = ~y,
        type = 'scatter', mode = 'markers',
        hovertext = ~tooltip, hoverinfo = 'text',
        marker = list(
          symbol = 'square',
          size = 7,
          color = '#333333',
          line = list(color = 'white', width = 1)
        ),
        name = 'Institutionen',
        showlegend = FALSE,
        inherit = FALSE
      )
    }

    # Field labels are placed near their clusters and repelled jointly
    # with person and institution labels.
    field_anchors <- net$centers |>
      transmute(label = field, kind = 'Field',
                x = cx, y = cy + radius * 0.85)
    labels <- bind_rows(
      nodes |> transmute(label, kind = if_else(type, 'Institution', 'Person'), x, y),
      field_anchors
    )
    all_x <- c(nodes$x, net$hulls$x, labels$x)
    all_y <- c(nodes$y, net$hulls$y, labels$y)
    xr <- range(all_x, finite = TRUE)
    yr <- range(all_y, finite = TRUE)
    xpad <- max(diff(xr) * 0.18, 4)
    ypad <- max(diff(yr) * 0.18, 4)
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)
    labels <- repel_labels(labels, xlim, ylim)

    # Small leader lines connect displaced labels to their nodes.
    displaced <- labels |>
      filter(kind != 'Field',
             sqrt(((lx - x) / diff(xlim))^2 +
                  ((ly - y) / diff(ylim))^2) > 0.018)
    if (nrow(displaced) > 0) {
      p <- p |> add_trace(
        x = as.vector(rbind(displaced$x, displaced$lx, NA_real_)),
        y = as.vector(rbind(displaced$y, displaced$ly, NA_real_)),
        type = 'scatter', mode = 'lines',
        line = list(color = 'rgba(90,90,90,0.35)', width = 0.7),
        hoverinfo = 'skip', showlegend = FALSE, inherit = FALSE
      )
    }

    # Permanent labels, positioned independently of Plotly's node traces.
    annotations <- lapply(seq_len(nrow(labels)), function(i) {
      row <- labels[i, ]
      list(
        x = row$lx, y = row$ly, text = htmltools::htmlEscape(row$label),
        showarrow = FALSE, xanchor = 'center', yanchor = 'middle',
        font = list(
          size = if (row$kind == 'Field') 14 else if (row$kind == 'Institution') 11 else 10,
          color = if (row$kind == 'Field') '#243746' else '#30343A',
          family = if (row$kind == 'Field') 'Arial Black, Arial' else 'Arial'
        )
      )
    })

    p |> layout(
      annotations = annotations,
      xaxis = list(visible = FALSE, zeroline = FALSE, range = xlim),
      yaxis = list(visible = FALSE, zeroline = FALSE,
                   scaleanchor = 'x', scaleratio = 1, range = ylim),
      hovermode = 'closest',
      legend = list(orientation = 'h', x = 0, y = -0.08),
      margin = list(l = 35, r = 35, t = 45, b = 75),
      paper_bgcolor = 'white', plot_bgcolor = 'white'
    ) |> config(displaylogo = FALSE, responsive = TRUE)
  })

  output$data <- renderTable({ data() }, striped = TRUE, spacing = 'xs')
}

shinyApp(ui = ui, server = server)
