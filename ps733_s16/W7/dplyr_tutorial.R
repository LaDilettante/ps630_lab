library(nycflights13)
library(dplyr)

flights

glimpse(flights)


filter() (and slice())
arrange()
select() (and rename())
distinct()
mutate() (and transmute())
summarise()
sample_n() and sample_frac()



flights %>%
  filter(month == 1)

flights %>%
  filter(month == 1 | month == 2, day == 1)



flights %>%
  group_by(year, month, day) %>%
  select(arr_delay, dep_delay) %>%
  summarise(
    arr = mean(arr_delay, na.rm = TRUE),
    dep = mean(dep_delay, na.rm = TRUE)
  ) %>%
  filter(arr > 30 | dep > 30)

# https://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html