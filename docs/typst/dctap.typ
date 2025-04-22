#import "@local/huangjac:1.0.0":*

#show: huang.with(
  title: [DC-TAP-seq Computational Analysis Tools],
  subtitle: [ Work In-progress],
  author: "Jacob Huang",
  date: datetime.today(),
  report-style: true,
)

#toc

#pagebreak()
#include "src/ch0-preface.typ"

#pagebreak()
#set page(numbering: "1")
#counter(page).update(1)
#include "src/qpcr.typ"
