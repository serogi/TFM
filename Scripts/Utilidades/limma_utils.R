### Función para obtener anovas de 2 vías con/sin interceptación a partir de Limma


Anova2way <- function(g1, g2, expr_data){
    design <- model.matrix(~ g1 * g2)
    rownames(design) <- colnames(expr_data)
    corte <- 0.05
    fit <- lmFit(expr_data, design)
    res <- eBayes(fit)
    res$p.adj <- apply(res$p.value, 2, p.adjust, method = "BH")
    results <- decideTest(res, p.value = corte)
    tt <- topTable(res, coef = 1:3)
    return(tt)
}

Anova2way2 <- function(g1, g2, expr_data){
}
