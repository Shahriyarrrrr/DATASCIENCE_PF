
install.packages(c("rentrez", "httr", "xml2", "tidyverse", "tidytext",
                   "SnowballC", "stopwords", "tm", "uwot", "dbscan",
                   "ggplot2", "ggrepel", "cluster", "aricode"))


library(rentrez)
library(httr)
library(xml2)
library(tidyverse)
library(tidytext)
library(SnowballC)
library(stopwords)
library(tm)
library(uwot)
library(dbscan)
library(ggplot2)
library(ggrepel)
library(cluster)
library(aricode)     


dir.create("output",       showWarnings = FALSE)
dir.create("output/data",  showWarnings = FALSE)
dir.create("output/plots", showWarnings = FALSE)


save_plot <- function(plot, filename, width = 10, height = 7) {
  path <- file.path("output", "plots", filename)
  ggsave(path, plot = plot, width = width, height = height, dpi = 300)
  cat("Plot saved:", path, "\n")
}


save_data <- function(data, filename) {
  path <- file.path("output", "data", filename)
  write_csv(data, path)
  cat("Data saved:", path, "\n")
}




fetch_pubmed <- function(query, n = 10000) {
  search_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
    "?db=pubmed",
    "&term=", URLencode(query),
    "&retmax=", n,
    "&usehistory=y",
    "&retmode=xml"
  )
  
  search_resp <- GET(search_url)
  search_xml  <- read_xml(content(search_resp, "text"))
  
  web_env   <- xml_text(xml_find_first(search_xml, "//WebEnv"))
  query_key <- xml_text(xml_find_first(search_xml, "//QueryKey"))
  count     <- as.integer(xml_text(xml_find_first(search_xml, "//Count")))
  
  cat("Total hits found:", count, "\n")
  cat("WebEnv:", web_env, "\n")
  cat("QueryKey:", query_key, "\n")
  
  total_to_fetch <- min(n, count)
  cat("Fetching:", total_to_fetch, "records...\n")
  
  batch_size   <- 200
  records_list <- list()
  
  for (start in seq(0, total_to_fetch - 1, by = batch_size)) {
    cat("Fetching records", start + 1, "to",
        min(start + batch_size, total_to_fetch), "...\n")
    
    fetch_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
      "?db=pubmed",
      "&query_key=", query_key,
      "&WebEnv=", web_env,
      "&retstart=", start,
      "&retmax=", batch_size,
      "&rettype=xml",
      "&retmode=xml"
    )
    
    fetch_resp <- GET(fetch_url)
    records_list[[length(records_list) + 1]] <- content(fetch_resp, "text")
    Sys.sleep(0.4)
  }
  
  parse_batch <- function(xml_text) {
    xml      <- read_xml(xml_text)
    articles <- xml_find_all(xml, "//PubmedArticle")
    
    if (length(articles) == 0) return(tibble())
    
    map_dfr(articles, function(article) {
      pmid      <- xml_text(xml_find_first(article, ".//PMID"))
      abs_nodes <- xml_find_all(article, ".//AbstractText")
      abstract  <- paste(xml_text(abs_nodes), collapse = " ")
      tibble(id = pmid, abstract = abstract, domain = 0L)
    })
  }
  
  map_dfr(records_list, parse_batch) %>%
    filter(nchar(abstract) > 50) %>%
    distinct(id, .keep_all = TRUE)
}

fetch_arxiv <- function(query, n = 10000) {
  batch_size   <- 200
  records_list <- list()
  
  for (start in seq(0, n - 1, by = batch_size)) {
    cat("Fetching arXiv records", start + 1, "to",
        min(start + batch_size, n), "...\n")
    
    url  <- paste0(
      "http://export.arxiv.org/api/query?search_query=all:",
      URLencode(query),
      "&start=", start,
      "&max_results=", batch_size
    )
    resp <- GET(url)
    records_list[[length(records_list) + 1]] <- content(resp, "text")
    Sys.sleep(3)
  }
  
  parse_arxiv_batch <- function(xml_text) {
    xml     <- read_xml(xml_text)
    ns      <- c(a = "http://www.w3.org/2005/Atom")
    entries <- xml_find_all(xml, "//a:entry", ns)
    
    if (length(entries) == 0) return(tibble())
    
    map_dfr(entries, function(entry) {
      id       <- xml_text(xml_find_first(entry, "a:id", ns))
      abstract <- xml_text(xml_find_first(entry, "a:summary", ns))
      tibble(id = id, abstract = abstract, domain = 1L)
    })
  }
  
  map_dfr(records_list, parse_arxiv_batch) %>%
    filter(nchar(abstract) > 50) %>%
    distinct(id, .keep_all = TRUE)
}

bio_df <- fetch_pubmed(
  query = "machine learning[tw] OR deep learning[tw] OR neural network[tw]",
  n     = 10000
)
save_data(bio_df, "01_bio_raw.csv")

ai_df <- fetch_arxiv(
  query = "deep learning OR neural network OR transformer",
  n     = 10000
)
save_data(ai_df, "01_ai_raw.csv")




combined_df <- bind_rows(bio_df, ai_df) %>%
  mutate(
    domain   = factor(domain, levels = c(0L, 1L), labels = c("Bio", "AI")),
    abstract = str_squish(abstract),
    doc_id   = row_number()
  ) %>%
  filter(nchar(abstract) > 100) %>%
  distinct(abstract, .keep_all = TRUE)

cat("Total documents:", nrow(combined_df), "\n")
cat("Bio:", sum(combined_df$domain == "Bio"),
    "| AI:", sum(combined_df$domain == "AI"), "\n")

save_data(combined_df, "02_combined.csv")





custom_stopwords <- c(
  stopwords::stopwords("en"),
  "study", "result", "results", "method", "methods", "approach",
  "using", "used", "use", "also", "may", "one", "two", "three",
  "however", "based", "proposed", "show", "shown", "paper",
  "model", "models", "data", "dataset", "task", "tasks"
)

tokens_df <- combined_df %>%
  select(doc_id, domain, abstract) %>%
  unnest_tokens(word, abstract) %>%
  filter(
    !word %in% custom_stopwords,
    str_detect(word, "^[a-z]{3,}$")
  ) %>%
  mutate(word = wordStem(word, language = "en")) %>%
  count(doc_id, domain, word, name = "n")

save_data(tokens_df, "03_tokens.csv")






tfidf_df <- tokens_df %>%
  bind_tf_idf(word, doc_id, n)

top_terms <- tfidf_df %>%
  group_by(word) %>%
  summarise(total = sum(tf_idf), .groups = "drop") %>%
  slice_max(total, n = 3000) %>%
  pull(word)

tfidf_filtered <- tfidf_df %>%
  filter(word %in% top_terms)

save_data(tfidf_filtered, "04_tfidf.csv")

dtm <- tfidf_filtered %>%
  cast_dtm(doc_id, word, tf_idf)

mat_doc_ids <- as.integer(rownames(dtm))
mat         <- as.matrix(dtm)

doc_meta <- combined_df %>%
  filter(doc_id %in% mat_doc_ids) %>%
  arrange(match(doc_id, mat_doc_ids))





set.seed(42)
umap_result <- umap(
  mat,
  n_components = 2,
  n_neighbors  = 15,
  min_dist     = 0.1,
  metric       = "cosine",
  verbose      = TRUE
)

umap_df <- doc_meta %>%
  mutate(
    umap1 = umap_result[, 1],
    umap2 = umap_result[, 2]
  )

save_data(umap_df, "05_umap.csv")





hdb <- hdbscan(umap_result, minPts = 10)
umap_df$cluster_hdb <- factor(hdb$cluster)

cat("HDBSCAN clusters found:", length(unique(hdb$cluster)) - 1,
    "(excluding noise)\n")
cat("Noise points (cluster 0):", sum(hdb$cluster == 0), "\n")


sil_scores <- map_dbl(2:10, function(k) {
  km  <- kmeans(umap_result, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist(umap_result))
  mean(sil[, 3])
})

best_k <- which.max(sil_scores) + 1
cat("Best k by silhouette:", best_k, "\n")

set.seed(42)
km <- kmeans(umap_result, centers = best_k, nstart = 25)
umap_df$cluster_km <- factor(km$cluster)


umap_df <- umap_df %>%
  mutate(cluster = cluster_km)

save_data(umap_df, "06_umap_clustered.csv")





cluster_composition <- umap_df %>%
  group_by(cluster, domain) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(
    total = sum(n),
    pct   = round(100 * n / total, 1)
  ) %>%
  ungroup()

cluster_type <- cluster_composition %>%
  pivot_wider(
    names_from  = domain,
    values_from = pct,
    values_fill = 0
  ) %>%
  group_by(cluster) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(type = case_when(
    Bio >= 70 ~ "Bio-dominant",
    AI  >= 70 ~ "AI-dominant",
    TRUE      ~ "Shared"
  ))

save_data(cluster_composition, "07_cluster_composition.csv")
save_data(cluster_type,        "07_cluster_type.csv")

print(cluster_type)






cat("\n=== Cluster Evaluation ===\n")

true_labels <- as.integer(umap_df$domain) - 1L
pred_labels <- as.integer(umap_df$cluster)


sil_obj   <- silhouette(km$cluster, dist(umap_result))
sil_score <- mean(sil_obj[, 3])
cat("Silhouette Score :", round(sil_score, 4),
    "(range -1 to +1, higher is better)\n")


sil_per_cluster <- as_tibble(sil_obj[, 1:3]) %>%
  rename(cluster = cluster, neighbor = neighbor, sil_width = sil_width) %>%
  group_by(cluster) %>%
  summarise(
    n        = n(),
    mean_sil = round(mean(sil_width), 4),
    min_sil  = round(min(sil_width),  4),
    max_sil  = round(max(sil_width),  4),
    .groups  = "drop"
  )

cat("\nPer-cluster Silhouette Scores:\n")
print(sil_per_cluster)
save_data(sil_per_cluster, "08_silhouette_per_cluster.csv")


nmi_score <- NMI(true_labels, pred_labels)
cat("\nNMI Score        :", round(nmi_score, 4),
    "(range 0 to 1, higher is better)\n")

ari_score <- ARI(true_labels, pred_labels)
cat("ARI Score        :", round(ari_score, 4),
    "(range -1 to +1, higher is better)\n")


eval_summary <- tibble(
  metric      = c("Silhouette Score", "NMI", "ARI"),
  value       = c(round(sil_score, 4),
                  round(nmi_score, 4),
                  round(ari_score, 4)),
  range       = c("-1 to +1", "0 to 1", "-1 to +1"),
  optimal     = c("Higher is better",
                  "Higher is better",
                  "Higher is better"),
  interpretation = c(
    case_when(
      sil_score >= 0.7  ~ "Strong cluster structure",
      sil_score >= 0.5  ~ "Reasonable cluster structure",
      sil_score >= 0.25 ~ "Weak cluster structure",
      TRUE              ~ "No substantial structure"
    ),
    case_when(
      nmi_score >= 0.75 ~ "High domain alignment",
      nmi_score >= 0.5  ~ "Moderate domain alignment",
      nmi_score >= 0.25 ~ "Low domain alignment",
      TRUE              ~ "Minimal domain alignment"
    ),
    case_when(
      ari_score >= 0.8  ~ "Excellent agreement",
      ari_score >= 0.6  ~ "Good agreement",
      ari_score >= 0.4  ~ "Moderate agreement",
      ari_score >= 0.2  ~ "Fair agreement",
      TRUE              ~ "Poor agreement"
    )
  )
)

cat("\n=== Evaluation Summary ===\n")
print(eval_summary)
save_data(eval_summary, "08_evaluation_summary.csv")


p_eval <- eval_summary %>%
  mutate(metric = factor(metric,
                         levels = c("Silhouette Score", "NMI", "ARI"))) %>%
  ggplot(aes(x = metric, y = value, fill = metric)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = paste0(value, "\n(", interpretation, ")")),
            vjust   = -0.4,
            size    = 3.2,
            lineheight = 0.9) +
  scale_fill_manual(values = c(
    "Silhouette Score" = "#1D9E75",
    "NMI"              = "#7F77DD",
    "ARI"              = "#D85A30"
  )) +
  scale_y_continuous(limits = c(0, max(eval_summary$value) * 1.4)) +
  theme_minimal(base_size = 13) +
  labs(
    title    = "Cluster Evaluation Metrics",
    subtitle = "Silhouette Score | NMI | ARI — all higher is better",
    x        = NULL,
    y        = "Score"
  )

save_plot(p_eval, "08_evaluation_metrics.png", width = 9, height = 6)


sil_k_df <- tibble(k = 2:10, silhouette_score = sil_scores)
save_data(sil_k_df, "08_silhouette_scores_by_k.csv")

p_sil_k <- ggplot(sil_k_df, aes(x = k, y = silhouette_score)) +
  geom_line(color = "#7F77DD", linewidth = 1) +
  geom_point(color = "#7F77DD", size = 3) +
  geom_vline(xintercept = best_k, linetype = "dashed",
             color = "#D85A30", linewidth = 0.8) +
  annotate("text", x = best_k + 0.2, y = min(sil_scores),
           label = paste("Best k =", best_k),
           hjust = 0, color = "#D85A30", size = 4) +
  theme_minimal(base_size = 12) +
  labs(title = "Silhouette scores by number of clusters (k)",
       x = "Number of clusters (k)",
       y = "Mean silhouette score")

save_plot(p_sil_k, "08_silhouette_by_k.png")


p_sil_cluster <- ggplot(sil_per_cluster,
                        aes(x    = factor(cluster),
                            y    = mean_sil,
                            fill = mean_sil)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = sil_score, linetype = "dashed",
             color = "#D85A30", linewidth = 0.8) +
  annotate("text", x = 0.6, y = sil_score + 0.005,
           label = paste("Overall mean =", round(sil_score, 3)),
           hjust = 0, color = "#D85A30", size = 3.5) +
  scale_fill_gradient(low = "#B5D4F4", high = "#185FA5") +
  theme_minimal(base_size = 12) +
  labs(title = "Mean silhouette score per cluster",
       x = "Cluster", y = "Mean silhouette score")

save_plot(p_sil_cluster, "08_silhouette_per_cluster.png")






cluster_keywords <- tfidf_filtered %>%
  inner_join(umap_df %>% select(doc_id, cluster), by = "doc_id") %>%
  group_by(cluster, word) %>%
  summarise(mean_tfidf = mean(tf_idf), .groups = "drop") %>%
  group_by(cluster) %>%
  slice_max(mean_tfidf, n = 10) %>%
  summarise(keywords = paste(word, collapse = ", "))

save_data(cluster_keywords, "09_cluster_keywords.csv")
print(cluster_keywords)





knowledge_map <- cluster_type %>%
  left_join(cluster_keywords, by = "cluster") %>%
  select(cluster, total, Bio, AI, type, keywords) %>%
  arrange(desc(type == "Shared"), desc(total))

save_data(knowledge_map, "10_knowledge_map.csv")

cat("\n=== Knowledge Map ===\n")
print(knowledge_map, n = Inf)

shared <- knowledge_map %>% filter(type == "Shared")
cat("\nShared clusters (AI-Bio overlap):\n")
print(shared)






cluster_centers <- umap_df %>%
  group_by(cluster) %>%
  summarise(umap1 = mean(umap1), umap2 = mean(umap2), .groups = "drop") %>%
  left_join(
    tfidf_filtered %>%
      inner_join(umap_df %>% select(doc_id, cluster), by = "doc_id") %>%
      group_by(cluster, word) %>%
      summarise(mean_tfidf = mean(tf_idf), .groups = "drop") %>%
      group_by(cluster) %>%
      slice_max(mean_tfidf, n = 3) %>%
      summarise(label = paste(word, collapse = "\n")),
    by = "cluster"
  )

umap_plot_df <- umap_df %>%
  left_join(cluster_type %>% select(cluster, type), by = "cluster")

p_umap <- ggplot(umap_plot_df, aes(umap1, umap2, color = cluster)) +
  geom_point(aes(shape = type), alpha = 0.4, size = 0.9) +
  geom_label_repel(
    data               = cluster_centers,
    aes(label          = label),
    size               = 2.8,
    label.size         = 0.25,
    label.padding      = unit(0.15, "lines"),
    box.padding        = unit(0.5,  "lines"),
    point.padding      = unit(0.3,  "lines"),
    max.overlaps       = Inf,
    min.segment.length = 0,
    show.legend        = FALSE,
    seed               = 42
  ) +
  scale_shape_manual(
    values = c("Bio-dominant" = 16, "AI-dominant" = 17, "Shared" = 15)
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  labs(
    title    = "UMAP projection by cluster",
    subtitle = paste0(
      "Silhouette = ", round(sil_score, 3),
      "  |  NMI = ",   round(nmi_score, 3),
      "  |  ARI = ",   round(ari_score, 3)
    ),
    x     = "UMAP 1",
    y     = "UMAP 2",
    color = "Cluster",
    shape = "Cluster type"
  )

save_plot(p_umap, "11_umap_by_cluster.png", width = 14, height = 9)

cat("\n=== All outputs saved ===\n")
cat("Data files : output/data/\n")
cat("Plot files : output/plots/\n")