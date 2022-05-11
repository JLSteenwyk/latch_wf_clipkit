<p align="center">
    <a href="https://latch.bio/">
        <img src="images/clipkit_logo.jpg" width=37.5%>
    </a>
    <a href="https://latch.bio/">
        <img src="images/latch_logo.png" width=50%>
    </a>
    <a href="https://console.latch.bio/explore/59878/info">
    <br/>
    <span style="font-size:larger;">Click here to see the workflow!</span></br>
    </a>
    </br></br>
    <a href="https://github.com/hannahle/alignment/graphs/contributors" alt="Contributors">
        <img src="https://img.shields.io/github/contributors/jlsteenwyk/latch_wf_clipkit">
    </a>
    <a href="https://twitter.com/intent/follow?screen_name=jlsteenwyk" alt="Author Twitter">
        <img src="https://img.shields.io/twitter/follow/jlsteenwyk?style=social&logo=twitter"
            alt="follow on Twitter">
    </a>
</p>

</br>

# ClipKIT, the multiple sequence alignment trimming toolkit
## About
ClipKIT is a fast and flexible alignment trimming
tool that keeps phylogenetically informative sites
and removes those that display characteristics poor
phylogenetic signal.

<br />

If you found clipkit useful, please cite *ClipKIT:
a multiple sequence alignment trimming software for
accurate phylogenomic inference*. Steenwyk et al. 2020,
PLoS Biology. doi:
[10.1371/journal.pbio.3001007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001007).

<br />

## Modes
Herein, we describe the various trimming modes implemented in ClipKIT. If you are unsure which is appropriate for you, we recommend using the default smart-gap trimming mode.

<br />

ClipKIT can be run with eight different modes, which are specified with the -m/–mode argument. Default: ‘smart-gap’
<br />
- smart-gap: dynamic determination of gaps threshold
- gappy: trim all sites that are above a threshold of gappyness (default: 0.9)
- kpic: keep only parismony informative and constant sites
- kpic-smart-gap: a combination of kpic- and smart-gap-based trimming
- kpic-gappy: a combination of kpic- and gappy-based trimming
- kpi: keep only parsimony informative sites
- kpi-smart-gap: a combination of kpi- and smart-gap-based trimming
- kpi-gappy: a combination of kpi- and gappy-based trimming
