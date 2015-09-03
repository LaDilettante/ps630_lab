sum=0
for (i in 1:6){
  for (j in 1:6){
    sum=i+j+sum
  }
}

sum/36

is <- rep(1:6, 6)
js <- rep(1:6, each = 6)

sum(is + js) / 36