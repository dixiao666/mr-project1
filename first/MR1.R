# 加载包
library("devtools")
library("ieugwasr")
library(TwoSampleMR)
library(ieugwasr)
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiIyMDI0MDA1MTA2NjhAc3htdS5lZHUuY24iLCJpYXQiOjE3NDIzNjkxNTUsImV4cCI6MTc0MzU3ODc1NX0.bSXTv8s-IL49_PMbPtw9WbqBa3BBO0uy3lgwR9TbRex0j52W16TTINBQFsmPNdmCdvx0Zv0J9bqSUSpwqmnz6vdEwFA-9lvUtdjBldAktrCWvTlgeL4HBzhVai3nFr62b0AA2xlztcycfbQdhd87KXDXeFqtR8XqxVZuQVpQQTCpiDL6ie0jwRg9lSC4PrNXzfvlFzlWcV7UNedQBunEwCfRgehcuUsQiwtTps2nzxQfoJjFnGgD5mw30y-2Bwf0rD1t8EiMZ4YOd6albLGZ904l0csKVoIfDfPvKG0rkl7jDAKZjBShjaNnssCUm_FhaID_3apFIyKv7ePRkjA-vw")
# 使用token
normalizePath("~/.Renviron")
ieugwasr::get_opengwas_jwt()
api_status()
user()
ieugwasr::user()
# 下载暴露数据
exposure_data <- extract_instruments("ukb-b-18336")
head(exposure_data)[1:4, ]

outcome_data <- extract_outcome_data(
  snps = exposure_data$SNP,       # 指定要提取结局数据的 SNP 列表
  outcomes = 'ukb-b-16890',       # 指定 GWAS 数据集
  proxies = FALSE,                # 不使用代理 SNP
  maf_threshold = 0.01,           # 最小等位基因频率阈值
)
data <- harmonise_data(
  exposure_dat = exposure_data,  # 指定暴露数据，这是之前从 GWAS 数据集中提取的工具变量数据
  outcome_dat = outcome_data)      # 指定结局数据，这是从另一个 GWAS 数据集中提取的关联结局数据
colnames(data)
# 正式开始 MR 分析啦！
result <- mr(data)
result

# 水平多效性检验
mr_pleiotropy_test(data)
# 散点图
p1 <- mr_scatter_plot(result, data)
p1

result_single <- mr_singlesnp(data)
p2 <- mr_forest_plot(result_single)
p2



