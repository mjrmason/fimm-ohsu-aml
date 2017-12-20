library("compareDrugScreens")
library("glmnet")

## This table has cell line doubling times (DoublingTime) for cell lines in CellLineName.
data("gcsi.line.info")

## Expression data (DESeq VST values) are here:
data("gcsi.genomics")
expr.cols <- colnames(gcsi.genomics)[grepl(pattern="vsd.GeneID", colnames(gcsi.genomics))]
gcsi.expr <- gcsi.genomics[, expr.cols]

common.cell.lines <- intersect(rownames(gcsi.line.info), rownames(gcsi.expr))
gcsi.expr <- gcsi.expr[common.cell.lines, ]
gcsi.line.info <- gcsi.line.info[common.cell.lines, ]

flag <- !is.na(gcsi.line.info$DoublingTime)
gcsi.line.info <- gcsi.line.info[flag, ]
gcsi.expr <- gcsi.expr[flag, ]

glm.fit <- cv.glmnet(x = gcsi.expr, y = gcsi.line.info$DoublingTime)

pred <- predict(glm.fit, newx = gcsi.expr)
x <- pred
y <- gcsi.line.info$DoublingTime
plot(x, y, xlim = c(min(x,y), max(x,y)), ylim = c(min(x,y), max(x,y)))
abline(a=0, b=1, lty=2)
cor.test(x, y)

lmf <- lm(y ~ x)
abline(lmf)

set.seed(1234)

common.cell.lines <- rownames(gcsi.expr)
train.index <- sample.int(n = length(common.cell.lines), size = floor(0.6 * length(common.cell.lines)))
train.cell.lines <- common.cell.lines[train.index]
test.cell.lines <- common.cell.lines[!(common.cell.lines %in% train.cell.lines)]

## alpha = 0: ridge
glm.fit <- cv.glmnet(x = gcsi.expr[train.cell.lines,], y = gcsi.line.info[train.cell.lines, "DoublingTime"], alpha = 1)

pred <- predict(glm.fit, newx = gcsi.expr[test.cell.lines,])
plot(pred, gcsi.line.info[test.cell.lines, "DoublingTime"])

flag <- pred < 100
actual <- gcsi.line.info[test.cell.lines, "DoublingTime"]
x <- pred[flag]
y <- actual[flag]
plot(x, y, xlim = c(min(x,y), max(x,y)), ylim = c(min(x,y), max(x,y)))
abline(a=0, b=1, lty=2)
cor.test(x, y)

lmf <- lm(y ~ x)
abline(lmf)



common.cell.lines <- rownames(gcsi.expr)
heme.line.info <- subset(gcsi.line.info, TissueMetaclass == "Blood")
common.cell.lines <- intersect(common.cell.lines, rownames(heme.line.info))
heme.expr <- gcsi.expr[common.cell.lines, ]
heme.line.info <- heme.line.info[common.cell.lines, ]
df <- as.data.frame(heme.expr)
df$SeedingDensity <- heme.line.info$SeedingDensity
heme.expr <- as.matrix(df)
train.index <- sample.int(n = length(common.cell.lines), size = floor(0.6 * length(common.cell.lines)))
train.cell.lines <- common.cell.lines[train.index]
test.cell.lines <- common.cell.lines[!(common.cell.lines %in% train.cell.lines)]

## alpha = 0: ridge
glm.fit <- cv.glmnet(x = heme.expr[train.cell.lines,], y = heme.line.info[train.cell.lines, "DoublingTime"])

pred <- predict(glm.fit, newx = heme.expr[test.cell.lines,])
plot(pred, heme.line.info[test.cell.lines, "DoublingTime"])

flag <- pred < 100
actual <- heme.line.info[test.cell.lines, "DoublingTime"]
x <- pred[flag]
y <- actual[flag]
smoothScatter(x, y, xlim = c(min(x,y), max(x,y)), ylim = c(min(x,y), max(x,y)))
abline(a=0, b=1, lty=2)
cor.test(x, y, method = "spearman")

lmf <- lm(y ~ x)
abline(lmf)

