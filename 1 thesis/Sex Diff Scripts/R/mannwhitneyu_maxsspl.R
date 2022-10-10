library(tidyverse)
library(ggplot2)

maxsspl = read.csv('D:/Tamara/Diskonnektionen/SSPLs/tidymaxSSPLs.csv')

wilcox.test(MaxSSPL~Sex, data = maxsspl, exact = FALSE, correct = FALSE, conf.int = FALSE)

  filter(Sex == "F")

maxsspl %>%
  group_by(Sex) %>%
  ggplot(
    aes(
      x = MaxSSPL,
      color = Sex,
      fill = Sex
    )
  ) + geom_bar(position = "dodge")
