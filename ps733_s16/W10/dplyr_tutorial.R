library(dplyr)
library(nycflights13)
data(flights)

dim(flights)

# Smart printing
flights

# Take a quick look
glimpse(flights)

# ---- 4 verbs ----
flights %>% select(year, month, day)

flights %>% filter(month == 2)

flights %>% arrange(dep_delay)
flights %>% arrange(desc(dep_delay))

flights %>% mutate(yearmonthday = paste(year, month, day, sep = "-")) %>%
  select(yearmonthday, year, month, day) %>%
  filter(month == 2) %>% summary()

filter(select(mutate()))

# ---- Annoying tasks that's really easy in dplyr ----

# Rename
# colnames(flights)[colnames(flights) == "dep_time"] == "Departure_Time"
flights %>%
  select(Depature_Time = dep_time)

# Select columns based on names
flights %>% select(starts_with("arr"))
flights %>% select(ends_with("time"))
flights %>% select(matches("rr"))

# Re-arrange columns
flights[ , c("dep_time", "dep_delay", "year", "month", "day")]
flights %>% select(dep_time, dep_delay, year, month, day, everything())

# group_by and summarise

# Chaining