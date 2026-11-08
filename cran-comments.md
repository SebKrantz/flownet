## Submission

This is a minor release (0.3.0 -> 0.3.1). All-or-nothing assignment now computes its shortest paths
in C, multithreaded with OpenMP; directed graph consolidation now also contracts two-way links given
as reciprocal arcs; and several bugs in network processing and parallel assignment are fixed. See
NEWS.md for details.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  Found the following (possibly) invalid URLs:
    URL: https://deepwiki.com/SebKrantz/flownet
      From: README.md
      Status: 429
      Message: Too Many Requests

  This is the target of the "Ask DeepWiki" badge in the README. The site rate-limits automated
  requests (HTTP 429); the URL is valid and opens normally in a browser. The same badge was present
  in flownet 0.3.0.

## Test environments

* local: macOS Sequoia 15.7 (aarch64-apple-darwin23), R 4.6.0, `R CMD check --as-cran`
* GitHub Actions: macOS (R release), Windows (R release), Ubuntu (R release and oldrel-1)

## Reverse dependencies

There are currently no reverse dependencies of flownet on CRAN.
