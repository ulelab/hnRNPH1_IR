awk -F'\t' -v tol=5 'BEGIN{OFS=FS}
NR==1 {print $0,"junction_type"; next}
{
  start_match = (($9-$4)<=tol && ($9-$4)>=-tol)
  end_match   = (($10-$5)<=tol && ($10-$5)>=-tol)

  jt="none"
  if (start_match && end_match) {
    jt="both"
  } else if (start_match) {
    jt = ($6=="+") ? "5ss" : "3ss"
  } else if (end_match) {
    jt = ($6=="+") ? "3ss" : "5ss"
  }

  print $0,jt
}' results/srajunctions.tsv > results/srajunctions.annotated.tsv
