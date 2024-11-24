create_heatmap <- function(df1, df2, top_n = 10, significance_threshold = 0.05) {
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(pheatmap)
  library(RColorBrewer)
  
  # GO 용어에서 숫자 제거 함수
  remove_go <- function(term) {
    str_remove(term, "\\s*\\(GO:\\d+\\)")
  }
  
  # 데이터 전처리
  df1 <- df1 %>%
    mutate(
      Condition = "R",
      logP = -log10(`P-value`),
      Term = remove_go(Term)
    )
  
  df2 <- df2 %>%
    mutate(
      Condition = "E",
      logP = -log10(`P-value`),
      Term = remove_go(Term)
    )
  
  # 상위 Term 선택
  top_terms_df1 <- df1 %>%
    arrange(`P-value`) %>%
    slice_head(n = top_n)
  
  top_terms_df2 <- df2 %>%
    arrange(`P-value`) %>%
    slice_head(n = top_n)
  
  selected_terms <- bind_rows(top_terms_df1, top_terms_df2) %>%
    distinct(Term)
  
  # 두 데이터프레임 병합
  combined_df <- bind_rows(df1, df2)
  
  # 선택된 Term 필터링
  combined_selected_df <- combined_df %>%
    filter(Term %in% selected_terms$Term)
  
  # 히트맵 데이터 준비 (wide format)
  heatmap_data <- combined_selected_df %>%
    select(Term, Condition, logP) %>%
    pivot_wider(
      names_from = Condition,
      values_from = logP,
      values_fill = list(logP = 0)
    )
  
  # P-value 데이터 준비 (wide format)
  pvalue_data <- combined_selected_df %>%
    select(Term, Condition, `P-value`) %>%
    pivot_wider(
      names_from = Condition,
      values_from = `P-value`,
      values_fill = list(`P-value` = 1)
    )
  
  # 행렬로 변환
  heatmap_matrix <- heatmap_data %>%
    column_to_rownames("Term") %>%
    as.matrix()
  
  pvalue_matrix <- pvalue_data %>%
    column_to_rownames("Term") %>%
    as.matrix()
  
  # 유의미성 표시 매트릭스 생성
  display_numbers <- ifelse(pvalue_matrix < significance_threshold, "*", "")
  
  # 히트맵 그리기
  pheatmap(
    heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    display_numbers = display_numbers,
    number_color = "black",
    color = brewer.pal(n = 9, name = "Reds"),
    main = "Heatmap of logP Values (Top Terms per Condition)",
    fontsize_number = 20,
    na_col = "white",
    angle_col = 0
  )
}