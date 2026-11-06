install.packages(c(
  "shiny", "dplyr", "readr",
  "ggplot2", "ggrepel", "plotly", "igraph"
))

library(shiny)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(igraph)

sheet_url <- paste0(
  "https://docs.google.com/spreadsheets/d/",
  "1y3qQ8eUz6yw8-pYHxwp87Dqy6Y3_cIdc6ZmbXjWCNcA",
  "/export?format=csv&gid=0"
)

ui <- fluidPage(
  titlePanel("Forschungslandschaft Schweiz"),
  p("Forschende nach Forschungsfeld und institutioneller Zugehörigkeit"),
  plotlyOutput("network", height = "850px"),
  h4("Aktuelle Daten"),
  tableOutput("data")
)

server <- function(input, output, session) {
  
  data <- reactive({
    invalidateLater(10000, session)
    
    df <- read_csv(sheet_url, show_col_types = FALSE)
    names(df) <- tolower(trimws(names(df)))
    
    required <- c(
      "institution", "person",
      "schwerpunkt", "forschungsfeld"
    )
    
    validate(need(
      all(required %in% names(df)),
      paste("Missing columns:",
            paste(setdiff(required, names(df)), collapse = ", "))
    ))
    
    df |>
      mutate(
        across(all_of(required), as.character),
        across(all_of(required), trimws),
        person = if_else(
          is.na(person) | person == "",
          "NA", person
        ),
        institution = if_else(
          is.na(institution) | institution == "",
          "Unbekannt", institution
        ),
        schwerpunkt = if_else(
          is.na(schwerpunkt) | schwerpunkt == "",
          "Nicht angegeben", schwerpunkt
        )
      ) |>
      filter(!is.na(forschungsfeld), forschungsfeld != "") |>
      distinct()
  })
  
  output$network <- renderPlotly({
    
    df <- data()
    validate(need(nrow(df) > 0, "Keine Daten vorhanden"))
    
    # Research-field cluster centres
    fields <- sort(unique(df$forschungsfeld))
    n_fields <- length(fields)
    
    angles <- seq(0, 2*pi, length.out = n_fields + 1)[
      -(n_fields + 1)
    ]
    
    ring_radius <- max(12, n_fields * 3.0)
    
    centres <- data.frame(
      forschungsfeld = fields,
      cx = ring_radius * cos(angles),
      cy = ring_radius * sin(angles)
    )
    
    # Position persons inside their research-field bubbles
    df <- df |>
      left_join(centres, by = "forschungsfeld") |>
      group_by(forschungsfeld) |>
      arrange(institution, person, .by_group = TRUE) |>
      mutate(
        node_number = row_number(),
        
        # Six nodes per ring, increasing capacity outward
        ring = ceiling(
          (sqrt(1 + 8 * node_number / 6) - 1) / 2
        ),
        
        previous_nodes = 3 * (ring - 1) * ring,
        
        index_in_ring = node_number - previous_nodes,
        
        nodes_in_ring = 6 * ring,
        
        theta = 2 * pi * (index_in_ring - 1) /
          nodes_in_ring,
        
        # Increased spacing between nodes
        radius = ring * 1.15,
        
        x = cx + radius * cos(theta),
        y = cy + radius * sin(theta)
      ) |>
      ungroup() |>
      mutate(node_id = row_number())
    
    
    # Create minimum spanning tree for each institution
    edge_list <- list()
    
    for (inst in unique(df$institution)) {
      
      sub <- df |> filter(institution == inst)
      
      if (nrow(sub) < 2) next
      
      coords <- as.matrix(sub[, c("x", "y")])
      distances <- as.matrix(dist(coords))
      
      graph <- graph_from_adjacency_matrix(
        distances,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )
      
      tree <- mst(graph)
      edges <- as_data_frame(tree, what = "edges")
      
      if (nrow(edges) == 0) next
      
      edge_list[[inst]] <- data.frame(
        x = sub$x[as.integer(edges$from)],
        y = sub$y[as.integer(edges$from)],
        xend = sub$x[as.integer(edges$to)],
        yend = sub$y[as.integer(edges$to)],
        institution = inst
      )
    }
    
    edges_df <- if (length(edge_list) > 0) {
      bind_rows(edge_list)
    } else {
      data.frame(
        x = numeric(), y = numeric(),
        xend = numeric(), yend = numeric(),
        institution = character()
      )
    }
    
    # Research-field bubble sizes
    # Larger research-field bubbles
    bubbles <- df |>
      group_by(forschungsfeld, cx, cy) |>
      summarise(
        radius = max(
          4.5,
          max(sqrt((x - cx)^2 + (y - cy)^2)) + 2.0
        ),
        .groups = "drop"
      )
    # Create circular boundaries
    circle_data <- bind_rows(lapply(
      seq_len(nrow(bubbles)),
      function(i) {
        a <- seq(0, 2*pi, length.out = 150)
        data.frame(
          forschungsfeld = bubbles$forschungsfeld[i],
          x = bubbles$cx[i] + bubbles$radius[i]*cos(a),
          y = bubbles$cy[i] + bubbles$radius[i]*sin(a)
        )
      }
    ))
    
    df <- df |>
      mutate(
        tooltip = paste0(
          person,
          "<br>Institution: ", institution,
          "<br>Schwerpunkt: ", schwerpunkt,
          "<br>Forschungsfeld: ", forschungsfeld
        )
      )
    
    p <- ggplot() +
      
      # Research-field bubbles
      geom_polygon(
        data = circle_data,
        aes(x = x, y = y, group = forschungsfeld),
        fill = "#DCEAF2",
        colour = "#A7C4D6",
        alpha = 0.35,
        linewidth = 0.5
      ) +
      
      # Institution edges
      geom_segment(
        data = edges_df,
        aes(
          x = x, y = y,
          xend = xend, yend = yend,
          colour = institution
        ),
        linewidth = 0.5,
        alpha = 0.4
      ) +
      
      # Person nodes
      geom_point(
        data = df,
        aes(
          x = x, y = y,
          colour = institution,
          text = tooltip
        ),
        size = 2.2
      ) +
      
      # Repelled person labels
      geom_text_repel(
        data = df,
        aes(x = x, y = y, label = person),
        size = 2.7,
        max.overlaps = Inf,
        box.padding = 0.7,
        point.padding = 0.4,
        force = 2,
        max.iter = 20000,
        segment.alpha = 0.3,
        seed = 42
      ) +
      
      # Research-field labels
      geom_text_repel(
        data = centres,
        aes(
          x = cx,
          y = cy - 3.8,
          label = forschungsfeld
        ),
        size = 2.7,
        max.overlaps = Inf,
        box.padding = 1,
        point.padding = 1,
        force = 4,
        max.iter = 20000,
        segment.alpha = 0.3,
        seed = 42
        fontface = "bold",
        size = 3.5
      ) +
      
      coord_equal(clip = "off") +
      theme_void(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.margin = margin(30, 50, 30, 50)
      ) +
      labs(colour = "Institution")
    
    ggplotly(p, tooltip = "text") |>
      layout(
        legend = list(orientation = "h"),
        hovermode = "closest"
      )
  })
  
  output$data <- renderTable(data())
}

shinyApp(ui, server)
