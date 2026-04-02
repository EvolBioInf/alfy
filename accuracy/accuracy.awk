BEGIN {
  if (!e || !o) {
    m = "Usage: awk -f accuracy.awk "
    m = m "-v e=<exp.txt> -v o=<obs.txt>"
    print m
    m = "Calculate accuracy from expected "
    m = m "and observed nearest neighbors."
    print m
    exit 0
  }
  files[0] = e
  files[1] = o
  for (i in files) {
    file = files[i]
    out = file ".bed"
    printf "" > out
    cmd = sprintf("cat %s", file)
    while (cmd | getline) {
      if ($1 ~ /^>/)
        continue
      start = $1 - 1
      end = $2
      sbjct = $4
      printf "%s\t%s\t%s\n",
        sbjct, start, end >> out
    }
    close(cmd)
    close(out)
  }
  e = files[0] ".bed"
  o = files[1] ".bed"
  cmd = "bedtools intersect -a %s -b %s"
  cmd = sprintf(cmd, e, o)
  while (cmd | getline) {
    l = $3 - $2
    tp += l
  }
  close(cmd)
  cmd = "cat %s"
  cmd = sprintf(cmd, e)
  c = 0
  while (cmd | getline) {
    c++
    if (c == 1)
      start = $2
    end = $3
  }
  close(cmd)
  len = end - start
  acc = tp / len
  printf("%.4f\n", acc)
  system("rm " e)
  system("rm " o)

}
