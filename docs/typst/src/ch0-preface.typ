#import "@local/huangjac:1.0.0":*

// Initial Preface Styles
#set heading(numbering: none)
#show heading: it => {
  block([
    #h(0.0em)
    #text(fill:colors.headers, it.body)
    #v(0.4em)
  ])
}
#show heading.where(level: 1): set text(size: 1.4em)

= Preface

== Purpose

Todo.

== Structure

Todo.

== Acknowledgements

Todo.
