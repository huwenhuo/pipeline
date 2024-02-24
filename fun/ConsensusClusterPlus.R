





<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
  <link rel="dns-prefetch" href="https://github.githubassets.com">
  <link rel="dns-prefetch" href="https://avatars0.githubusercontent.com">
  <link rel="dns-prefetch" href="https://avatars1.githubusercontent.com">
  <link rel="dns-prefetch" href="https://avatars2.githubusercontent.com">
  <link rel="dns-prefetch" href="https://avatars3.githubusercontent.com">
  <link rel="dns-prefetch" href="https://github-cloud.s3.amazonaws.com">
  <link rel="dns-prefetch" href="https://user-images.githubusercontent.com/">



  <link crossorigin="anonymous" media="all" integrity="sha512-G5IW3QX9jLeIufJaob0LkBXhXTZSiWUoXGNTvU9KgK4dfhMjKB3hfAy9hjsk5YYoN0GA3b0yekVqlMA5EYlDDA==" rel="stylesheet" href="https://github.githubassets.com/assets/frameworks-c567bfeb1cb9f4ac89533a5e03dbe623.css" />
  <link crossorigin="anonymous" media="all" integrity="sha512-5LDBQ/DhsFq/ry7aepdbc1FV0mW/0lPqH6GlYFfzzd5wo8MBz9Rb67LigKe0uKM8kK5QrA/QsJCpsE1+hOLF6w==" rel="stylesheet" href="https://github.githubassets.com/assets/site-0e25450ff46c1fbf61d504f9ffc4d784.css" />
    <link crossorigin="anonymous" media="all" integrity="sha512-5v4G7h0D7l4LzJAV5ZDxzCmB2bqr6u6GwZOBuLnpxLeh5tzwgG4JHHfIW05l3333JLycXT05cHmTgkwC9I4fNQ==" rel="stylesheet" href="https://github.githubassets.com/assets/github-89642c25e1c39ca29fd1d22a91c65558.css" />
    
    
    
    

  <meta name="viewport" content="width=device-width">
  
  <title>ConsensusClusterPlus/ConsensusClusterPlus.R at master · renzhonglu/ConsensusClusterPlus · GitHub</title>
    <meta name="description" content="This is a read-only mirror of the Bioconductor SVN repository. Package Homepage: http://bioconductor.org/packages/devel/bioc/html/ConsensusClusterPlus.html Bug Reports: https://support.bioconductor.org/p/new/post/?tag_val=ConsensusClusterPlus. - renzhonglu/ConsensusClusterPlus">
    <link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="GitHub">
  <link rel="fluid-icon" href="https://github.com/fluidicon.png" title="GitHub">
  <meta property="fb:app_id" content="1401488693436528">

    
    <meta property="og:image" content="https://avatars3.githubusercontent.com/u/8721120?s=400&amp;v=4" /><meta property="og:site_name" content="GitHub" /><meta property="og:type" content="object" /><meta property="og:title" content="renzhonglu/ConsensusClusterPlus" /><meta property="og:url" content="https://github.com/renzhonglu/ConsensusClusterPlus" /><meta property="og:description" content="This is a read-only mirror of the Bioconductor SVN repository. Package Homepage: http://bioconductor.org/packages/devel/bioc/html/ConsensusClusterPlus.html Bug Reports: https://support.bioconductor..." />

  <link rel="assets" href="https://github.githubassets.com/">
  
  <meta name="pjax-timeout" content="1000">
  
  <meta name="request-id" content="CFC9:3264:93143:10807E:5C911135" data-pjax-transient>


  

  <meta name="selected-link" value="repo_source" data-pjax-transient>

      <meta name="google-site-verification" content="KT5gs8h0wvaagLKAVWq8bbeNwnZZK1r1XQysX3xurLU">
    <meta name="google-site-verification" content="ZzhVyEFwb7w3e0-uOTltm8Jsck2F5StVihD0exw2fsA">
    <meta name="google-site-verification" content="GXs5KoUUkNCoaAZn7wPN-t01Pywp9M3sEjnt_3_ZWPc">

  <meta name="octolytics-host" content="collector.githubapp.com" /><meta name="octolytics-app-id" content="github" /><meta name="octolytics-event-url" content="https://collector.githubapp.com/github-external/browser_event" /><meta name="octolytics-dimension-request_id" content="CFC9:3264:93143:10807E:5C911135" /><meta name="octolytics-dimension-region_edge" content="iad" /><meta name="octolytics-dimension-region_render" content="iad" />
<meta name="analytics-location" content="/&lt;user-name&gt;/&lt;repo-name&gt;/blob/show" data-pjax-transient="true" />



    <meta name="google-analytics" content="UA-3769691-2">


<meta class="js-ga-set" name="dimension1" content="Logged Out">



  

      <meta name="hostname" content="github.com">
    <meta name="user-login" content="">

      <meta name="expected-hostname" content="github.com">
    <meta name="js-proxy-site-detection-payload" content="MTMwY2ZlODM1ZTRjZTEyMmQ1M2NmZjQ1Mzc1ODNkYzc1ZjRlMzI3ZGQzMWYyZTIxNWJhZmNkZTA5NDJlZWI1ZHx7InJlbW90ZV9hZGRyZXNzIjoiMTQwLjE2My43MC44IiwicmVxdWVzdF9pZCI6IkNGQzk6MzI2NDo5MzE0MzoxMDgwN0U6NUM5MTExMzUiLCJ0aW1lc3RhbXAiOjE1NTMwMTA5OTgsImhvc3QiOiJnaXRodWIuY29tIn0=">

    <meta name="enabled-features" content="UNIVERSE_BANNER,MARKETPLACE_SOCIAL_PROOF,MARKETPLACE_PLAN_RESTRICTION_EDITOR">

  <meta name="html-safe-nonce" content="3084c947833c4c02b153a614e19e488a5223f0fb">

  <meta http-equiv="x-pjax-version" content="654cdaf3b791c3b18a14dc0a447dceb4">
  

      <link href="https://github.com/renzhonglu/ConsensusClusterPlus/commits/master.atom" rel="alternate" title="Recent Commits to ConsensusClusterPlus:master" type="application/atom+xml">

  <meta name="go-import" content="github.com/renzhonglu/ConsensusClusterPlus git https://github.com/renzhonglu/ConsensusClusterPlus.git">

  <meta name="octolytics-dimension-user_id" content="8721120" /><meta name="octolytics-dimension-user_login" content="renzhonglu" /><meta name="octolytics-dimension-repository_id" content="39429778" /><meta name="octolytics-dimension-repository_nwo" content="renzhonglu/ConsensusClusterPlus" /><meta name="octolytics-dimension-repository_public" content="true" /><meta name="octolytics-dimension-repository_is_fork" content="false" /><meta name="octolytics-dimension-repository_network_root_id" content="39429778" /><meta name="octolytics-dimension-repository_network_root_nwo" content="renzhonglu/ConsensusClusterPlus" /><meta name="octolytics-dimension-repository_explore_github_marketplace_ci_cta_shown" content="false" />


    <link rel="canonical" href="https://github.com/renzhonglu/ConsensusClusterPlus/blob/master/R/ConsensusClusterPlus.R" data-pjax-transient>


  <meta name="browser-stats-url" content="https://api.github.com/_private/browser/stats">

  <meta name="browser-errors-url" content="https://api.github.com/_private/browser/errors">

  <link rel="mask-icon" href="https://github.githubassets.com/pinned-octocat.svg" color="#000000">
  <link rel="icon" type="image/x-icon" class="js-site-favicon" href="https://github.githubassets.com/favicon.ico">

<meta name="theme-color" content="#1e2327">




  <link rel="manifest" href="/manifest.json" crossOrigin="use-credentials">

  </head>

  <body class="logged-out env-production page-blob">
    

  <div class="position-relative js-header-wrapper ">
    <a href="#start-of-content" tabindex="1" class="px-2 py-4 bg-blue text-white show-on-focus js-skip-to-content">Skip to content</a>
    <div id="js-pjax-loader-bar" class="pjax-loader-bar"><div class="progress"></div></div>

    
    
    


        
<header class="Header header-logged-out  position-relative f4 py-3" role="banner">
  <div class="container-lg d-flex px-3">
    <div class="d-flex flex-justify-between flex-items-center">
        <a class="mr-4" href="https://github.com/" aria-label="Homepage" data-ga-click="(Logged out) Header, go to homepage, icon:logo-wordmark">
          <svg height="32" class="octicon octicon-mark-github text-white" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
        </a>
    </div>

    <div class="HeaderMenu HeaderMenu--logged-out d-flex flex-justify-between flex-items-center flex-auto">
      <div class="d-none">
        <button class="btn-link js-details-target" type="button" aria-label="Toggle navigation" aria-expanded="false">
          <svg height="24" class="octicon octicon-x text-gray" viewBox="0 0 12 16" version="1.1" width="18" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
        </button>
      </div>

        <nav class="mt-0" aria-label="Global">
          <ul class="d-flex list-style-none">
              <li class=" mr-3 mr-lg-3 edge-item-fix position-relative flex-wrap flex-justify-between d-flex flex-items-center ">
                <details class="HeaderMenu-details details-overlay details-reset width-full">
                  <summary class="HeaderMenu-summary HeaderMenu-link px-0 py-3 border-0 no-wrap  d-inline-block">
                    Why GitHub?
                    <svg x="0px" y="0px" viewBox="0 0 14 8" xml:space="preserve" fill="none" class="icon-chevon-down-mktg position-relative">
                      <path d="M1,1l6.2,6L13,1"></path>
                    </svg>
                  </summary>
                  <div class="dropdown-menu flex-auto rounded-1 bg-white px-0 mt-0  p-4 left-n4 position-absolute">
                    <a href="/features" class="py-2 lh-condensed-ultra d-block link-gray-dark no-underline h5 Bump-link--hover" data-ga-click="(Logged out) Header, go to Features">Features <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a>
                    <ul class="list-style-none f5 pb-3">
                      <li class="edge-item-fix"><a href="/features/code-review/" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Code review">Code review</a></li>
                      <li class="edge-item-fix"><a href="/features/project-management/" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Project management">Project management</a></li>
                      <li class="edge-item-fix"><a href="/features/integrations" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Integrations">Integrations</a></li>
                      <li class="edge-item-fix"><a href="/features/actions" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Actions">Actions</a>
                      <li class="edge-item-fix"><a href="/features#team-management" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Team management">Team management</a></li>
                      <li class="edge-item-fix"><a href="/features#social-coding" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Social coding">Social coding</a></li>
                      <li class="edge-item-fix"><a href="/features#documentation" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Documentation">Documentation</a></li>
                      <li class="edge-item-fix"><a href="/features#code-hosting" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Code hosting">Code hosting</a></li>
                    </ul>

                    <ul class="list-style-none mb-0 border-lg-top pt-lg-3">
                      <li class="edge-item-fix"><a href="/case-studies" class="py-2 lh-condensed-ultra d-block no-underline link-gray-dark no-underline h5 Bump-link--hover" data-ga-click="(Logged out) Header, go to Customer stories">Customer stories <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a></li>
                      <li class="edge-item-fix"><a href="/security" class="py-2 lh-condensed-ultra d-block no-underline link-gray-dark no-underline h5 Bump-link--hover" data-ga-click="(Logged out) Header, go to Security">Security <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a></li>
                    </ul>
                  </div>
                </details>
              </li>
              <li class=" mr-3 mr-lg-3">
                <a href="/enterprise" class="HeaderMenu-link no-underline py-3 d-block d-lg-inline-block" data-ga-click="(Logged out) Header, click, go to Enterprise">Enterprise</a>
              </li>

              <li class=" mr-3 mr-lg-3 edge-item-fix position-relative flex-wrap flex-justify-between d-flex flex-items-center ">
                <details class="HeaderMenu-details details-overlay details-reset width-full">
                  <summary class="HeaderMenu-summary HeaderMenu-link px-0 py-3 border-0 no-wrap  d-inline-block">
                    Explore
                    <svg x="0px" y="0px" viewBox="0 0 14 8" xml:space="preserve" fill="none" class="icon-chevon-down-mktg position-relative">
                      <path d="M1,1l6.2,6L13,1"></path>
                    </svg>
                  </summary>

                  <div class="dropdown-menu flex-auto rounded-1 bg-white px-0 pt-2 pb-0 mt-0  p-4 left-n4 position-absolute">
                    <ul class="list-style-none mb-3">
                      <li class="edge-item-fix"><a href="/explore" class="py-2 lh-condensed-ultra d-block link-gray-dark no-underline h5 Bump-link--hover" data-ga-click="(Logged out) Header, go to Features">Explore GitHub <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a></li>
                    </ul>

                    <h4 class="text-gray-light text-normal text-mono f5 mb-2  border-top pt-3">Learn &amp; contribute</h4>
                    <ul class="list-style-none mb-3">
                      <li class="edge-item-fix"><a href="/topics" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Topics">Topics</a></li>
                        <li class="edge-item-fix"><a href="/collections" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Collections">Collections</a></li>
                      <li class="edge-item-fix"><a href="/trending" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Trending">Trending</a></li>
                      <li class="edge-item-fix"><a href="https://lab.github.com/" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Learning lab">Learning Lab</a></li>
                      <li class="edge-item-fix"><a href="https://opensource.guide" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Open source guides">Open source guides</a></li>
                    </ul>

                    <h4 class="text-gray-light text-normal text-mono f5 mb-2  border-top pt-3">Connect with others</h4>
                    <ul class="list-style-none mb-0">
                      <li class="edge-item-fix"><a href="https://github.com/events" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Events">Events</a></li>
                      <li class="edge-item-fix"><a href="https://github.community" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Community forum">Community forum</a></li>
                      <li class="edge-item-fix"><a href="https://education.github.com" class="py-2 pb-0 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to GitHub Education">GitHub Education</a></li>
                    </ul>
                  </div>
                </details>
              </li>

              <li class=" mr-3 mr-lg-3">
                <a href="/marketplace" class="HeaderMenu-link no-underline py-3 d-block d-lg-inline-block" data-ga-click="(Logged out) Header, go to Marketplace">Marketplace</a>
              </li>

              <li class=" mr-3 mr-lg-3 edge-item-fix position-relative flex-wrap flex-justify-between d-flex flex-items-center ">
                <details class="HeaderMenu-details details-overlay details-reset width-full">
                  <summary class="HeaderMenu-summary HeaderMenu-link px-0 py-3 border-0 no-wrap  d-inline-block">
                    Pricing
                    <svg x="0px" y="0px" viewBox="0 0 14 8" xml:space="preserve" fill="none" class="icon-chevon-down-mktg position-relative">
                       <path d="M1,1l6.2,6L13,1"></path>
                    </svg>
                  </summary>

                  <div class="dropdown-menu flex-auto rounded-1 bg-white px-0 pt-2 pb-4 mt-0  p-4 left-n4 position-absolute">
                    <a href="/pricing" class="pb-2 lh-condensed-ultra d-block link-gray-dark no-underline h5 Bump-link--hover" data-ga-click="(Logged out) Header, go to Pricing">Plans <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a>

                    <ul class="list-style-none mb-3">
                      <li class="edge-item-fix"><a href="/pricing#feature-comparison" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Compare features">Compare plans</a></li>
                      <li class="edge-item-fix"><a href="https://enterprise.github.com/contact" class="py-2 lh-condensed-ultra d-block link-gray no-underline f5" data-ga-click="(Logged out) Header, go to Compare features">Contact Sales</a></li>
                    </ul>

                    <ul class="list-style-none mb-0  border-top pt-3">
                      <li class="edge-item-fix"><a href="/nonprofit" class="py-2 lh-condensed-ultra d-block no-underline link-gray-dark no-underline h5 Bump-link--hover" data-ga-click="(Logged out) Header, go to Nonprofits">Nonprofit <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a></li>
                      <li class="edge-item-fix"><a href="https://education.github.com" class="py-2 pb-0 lh-condensed-ultra d-block no-underline link-gray-dark no-underline h5 Bump-link--hover"  data-ga-click="(Logged out) Header, go to Education">Education <span class="Bump-link-symbol float-right text-normal text-gray-light">&rarr;</span></a></li>
                    </ul>
                  </div>
                </details>
              </li>
          </ul>
        </nav>

      <div class="d-flex flex-items-center px-0 text-center text-left">
          <div class="d-lg-flex mr-3">
            <div class="header-search scoped-search site-scoped-search js-site-search position-relative js-jump-to"
  role="combobox"
  aria-owns="jump-to-results"
  aria-label="Search or jump to"
  aria-haspopup="listbox"
  aria-expanded="false"
>
  <div class="position-relative">
    <!-- '"` --><!-- </textarea></xmp> --></option></form><form class="js-site-search-form" role="search" aria-label="Site" data-scope-type="Repository" data-scope-id="39429778" data-scoped-search-url="/renzhonglu/ConsensusClusterPlus/search" data-unscoped-search-url="/search" action="/renzhonglu/ConsensusClusterPlus/search" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" />
      <label class="form-control header-search-wrapper header-search-wrapper-jump-to position-relative d-flex flex-justify-between flex-items-center js-chromeless-input-container">
        <input type="text"
          class="form-control header-search-input jump-to-field js-jump-to-field js-site-search-focus js-site-search-field is-clearable"
          data-hotkey="s,/"
          name="q"
          value=""
          placeholder="Search"
          data-unscoped-placeholder="Search GitHub"
          data-scoped-placeholder="Search"
          autocapitalize="off"
          aria-autocomplete="list"
          aria-controls="jump-to-results"
          aria-label="Search"
          data-jump-to-suggestions-path="/_graphql/GetSuggestedNavigationDestinations#csrf-token=NIwip6z+CUhn/wbemVqdQ4YudUF23NuLVttTCXiPWZ657BYtgKS+0FWseE9K/YHZh+1C/VZ4bJqKXO3NveoKLw=="
          spellcheck="false"
          autocomplete="off"
          >
          <input type="hidden" class="js-site-search-type-field" name="type" >
            <img src="https://github.githubassets.com/images/search-key-slash.svg" alt="" class="mr-2 header-search-key-slash">

            <div class="Box position-absolute overflow-hidden d-none jump-to-suggestions js-jump-to-suggestions-container">
              
<ul class="d-none js-jump-to-suggestions-template-container">
  

<li class="d-flex flex-justify-start flex-items-center p-0 f5 navigation-item js-navigation-item js-jump-to-suggestion" role="option">
  <a tabindex="-1" class="no-underline d-flex flex-auto flex-items-center jump-to-suggestions-path js-jump-to-suggestion-path js-navigation-open p-2" href="">
    <div class="jump-to-octicon js-jump-to-octicon flex-shrink-0 mr-2 text-center d-none">
      <svg height="16" width="16" class="octicon octicon-repo flex-shrink-0 js-jump-to-octicon-repo d-none" title="Repository" aria-label="Repository" viewBox="0 0 12 16" version="1.1" role="img"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
      <svg height="16" width="16" class="octicon octicon-project flex-shrink-0 js-jump-to-octicon-project d-none" title="Project" aria-label="Project" viewBox="0 0 15 16" version="1.1" role="img"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
      <svg height="16" width="16" class="octicon octicon-search flex-shrink-0 js-jump-to-octicon-search d-none" title="Search" aria-label="Search" viewBox="0 0 16 16" version="1.1" role="img"><path fill-rule="evenodd" d="M15.7 13.3l-3.81-3.83A5.93 5.93 0 0 0 13 6c0-3.31-2.69-6-6-6S1 2.69 1 6s2.69 6 6 6c1.3 0 2.48-.41 3.47-1.11l3.83 3.81c.19.2.45.3.7.3.25 0 .52-.09.7-.3a.996.996 0 0 0 0-1.41v.01zM7 10.7c-2.59 0-4.7-2.11-4.7-4.7 0-2.59 2.11-4.7 4.7-4.7 2.59 0 4.7 2.11 4.7 4.7 0 2.59-2.11 4.7-4.7 4.7z"/></svg>
    </div>

    <img class="avatar mr-2 flex-shrink-0 js-jump-to-suggestion-avatar d-none" alt="" aria-label="Team" src="" width="28" height="28">

    <div class="jump-to-suggestion-name js-jump-to-suggestion-name flex-auto overflow-hidden text-left no-wrap css-truncate css-truncate-target">
    </div>

    <div class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none js-jump-to-badge-search">
      <span class="js-jump-to-badge-search-text-default d-none" aria-label="in this repository">
        In this repository
      </span>
      <span class="js-jump-to-badge-search-text-global d-none" aria-label="in all of GitHub">
        All GitHub
      </span>
      <span aria-hidden="true" class="d-inline-block ml-1 v-align-middle">↵</span>
    </div>

    <div aria-hidden="true" class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none d-on-nav-focus js-jump-to-badge-jump">
      Jump to
      <span class="d-inline-block ml-1 v-align-middle">↵</span>
    </div>
  </a>
</li>

</ul>

<ul class="d-none js-jump-to-no-results-template-container">
  <li class="d-flex flex-justify-center flex-items-center f5 d-none js-jump-to-suggestion p-2">
    <span class="text-gray">No suggested jump to results</span>
  </li>
</ul>

<ul id="jump-to-results" role="listbox" class="p-0 m-0 js-navigation-container jump-to-suggestions-results-container js-jump-to-suggestions-results-container">
  

<li class="d-flex flex-justify-start flex-items-center p-0 f5 navigation-item js-navigation-item js-jump-to-scoped-search d-none" role="option">
  <a tabindex="-1" class="no-underline d-flex flex-auto flex-items-center jump-to-suggestions-path js-jump-to-suggestion-path js-navigation-open p-2" href="">
    <div class="jump-to-octicon js-jump-to-octicon flex-shrink-0 mr-2 text-center d-none">
      <svg height="16" width="16" class="octicon octicon-repo flex-shrink-0 js-jump-to-octicon-repo d-none" title="Repository" aria-label="Repository" viewBox="0 0 12 16" version="1.1" role="img"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
      <svg height="16" width="16" class="octicon octicon-project flex-shrink-0 js-jump-to-octicon-project d-none" title="Project" aria-label="Project" viewBox="0 0 15 16" version="1.1" role="img"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
      <svg height="16" width="16" class="octicon octicon-search flex-shrink-0 js-jump-to-octicon-search d-none" title="Search" aria-label="Search" viewBox="0 0 16 16" version="1.1" role="img"><path fill-rule="evenodd" d="M15.7 13.3l-3.81-3.83A5.93 5.93 0 0 0 13 6c0-3.31-2.69-6-6-6S1 2.69 1 6s2.69 6 6 6c1.3 0 2.48-.41 3.47-1.11l3.83 3.81c.19.2.45.3.7.3.25 0 .52-.09.7-.3a.996.996 0 0 0 0-1.41v.01zM7 10.7c-2.59 0-4.7-2.11-4.7-4.7 0-2.59 2.11-4.7 4.7-4.7 2.59 0 4.7 2.11 4.7 4.7 0 2.59-2.11 4.7-4.7 4.7z"/></svg>
    </div>

    <img class="avatar mr-2 flex-shrink-0 js-jump-to-suggestion-avatar d-none" alt="" aria-label="Team" src="" width="28" height="28">

    <div class="jump-to-suggestion-name js-jump-to-suggestion-name flex-auto overflow-hidden text-left no-wrap css-truncate css-truncate-target">
    </div>

    <div class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none js-jump-to-badge-search">
      <span class="js-jump-to-badge-search-text-default d-none" aria-label="in this repository">
        In this repository
      </span>
      <span class="js-jump-to-badge-search-text-global d-none" aria-label="in all of GitHub">
        All GitHub
      </span>
      <span aria-hidden="true" class="d-inline-block ml-1 v-align-middle">↵</span>
    </div>

    <div aria-hidden="true" class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none d-on-nav-focus js-jump-to-badge-jump">
      Jump to
      <span class="d-inline-block ml-1 v-align-middle">↵</span>
    </div>
  </a>
</li>

  

<li class="d-flex flex-justify-start flex-items-center p-0 f5 navigation-item js-navigation-item js-jump-to-global-search d-none" role="option">
  <a tabindex="-1" class="no-underline d-flex flex-auto flex-items-center jump-to-suggestions-path js-jump-to-suggestion-path js-navigation-open p-2" href="">
    <div class="jump-to-octicon js-jump-to-octicon flex-shrink-0 mr-2 text-center d-none">
      <svg height="16" width="16" class="octicon octicon-repo flex-shrink-0 js-jump-to-octicon-repo d-none" title="Repository" aria-label="Repository" viewBox="0 0 12 16" version="1.1" role="img"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
      <svg height="16" width="16" class="octicon octicon-project flex-shrink-0 js-jump-to-octicon-project d-none" title="Project" aria-label="Project" viewBox="0 0 15 16" version="1.1" role="img"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
      <svg height="16" width="16" class="octicon octicon-search flex-shrink-0 js-jump-to-octicon-search d-none" title="Search" aria-label="Search" viewBox="0 0 16 16" version="1.1" role="img"><path fill-rule="evenodd" d="M15.7 13.3l-3.81-3.83A5.93 5.93 0 0 0 13 6c0-3.31-2.69-6-6-6S1 2.69 1 6s2.69 6 6 6c1.3 0 2.48-.41 3.47-1.11l3.83 3.81c.19.2.45.3.7.3.25 0 .52-.09.7-.3a.996.996 0 0 0 0-1.41v.01zM7 10.7c-2.59 0-4.7-2.11-4.7-4.7 0-2.59 2.11-4.7 4.7-4.7 2.59 0 4.7 2.11 4.7 4.7 0 2.59-2.11 4.7-4.7 4.7z"/></svg>
    </div>

    <img class="avatar mr-2 flex-shrink-0 js-jump-to-suggestion-avatar d-none" alt="" aria-label="Team" src="" width="28" height="28">

    <div class="jump-to-suggestion-name js-jump-to-suggestion-name flex-auto overflow-hidden text-left no-wrap css-truncate css-truncate-target">
    </div>

    <div class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none js-jump-to-badge-search">
      <span class="js-jump-to-badge-search-text-default d-none" aria-label="in this repository">
        In this repository
      </span>
      <span class="js-jump-to-badge-search-text-global d-none" aria-label="in all of GitHub">
        All GitHub
      </span>
      <span aria-hidden="true" class="d-inline-block ml-1 v-align-middle">↵</span>
    </div>

    <div aria-hidden="true" class="border rounded-1 flex-shrink-0 bg-gray px-1 text-gray-light ml-1 f6 d-none d-on-nav-focus js-jump-to-badge-jump">
      Jump to
      <span class="d-inline-block ml-1 v-align-middle">↵</span>
    </div>
  </a>
</li>


</ul>

            </div>
      </label>
</form>  </div>
</div>

          </div>

        <a class="HeaderMenu-link no-underline mr-3" href="/login?return_to=%2Frenzhonglu%2FConsensusClusterPlus%2Fblob%2Fmaster%2FR%2FConsensusClusterPlus.R" data-ga-click="(Logged out) Header, clicked Sign in, text:sign-in">Sign&nbsp;in</a>
          <a class="HeaderMenu-link d-inline-block no-underline border border-gray-dark rounded-1 px-2 py-1" href="/join" data-ga-click="(Logged out) Header, clicked Sign up, text:sign-up">Sign&nbsp;up</a>
      </div>
    </div>
  </div>
</header>

  </div>

  <div id="start-of-content" class="show-on-focus"></div>

    <div id="js-flash-container">

</div>



  <div class="application-main " data-commit-hovercards-enabled>
        <div itemscope itemtype="http://schema.org/SoftwareSourceCode" class="">
    <main id="js-repo-pjax-container" data-pjax-container >
      


  






  <div class="pagehead repohead instapaper_ignore readability-menu experiment-repo-nav  ">
    <div class="repohead-details-container clearfix container">

      <ul class="pagehead-actions">



  <li>
      <a href="/login?return_to=%2Frenzhonglu%2FConsensusClusterPlus"
    class="btn btn-sm btn-with-count tooltipped tooltipped-s"
    aria-label="You must be signed in to watch a repository" rel="nofollow">
    <svg class="octicon octicon-eye v-align-text-bottom" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8.06 2C3 2 0 8 0 8s3 6 8.06 6C13 14 16 8 16 8s-3-6-7.94-6zM8 12c-2.2 0-4-1.78-4-4 0-2.2 1.8-4 4-4 2.22 0 4 1.8 4 4 0 2.22-1.78 4-4 4zm2-4c0 1.11-.89 2-2 2-1.11 0-2-.89-2-2 0-1.11.89-2 2-2 1.11 0 2 .89 2 2z"/></svg>
    Watch
  </a>
  <a class="social-count" href="/renzhonglu/ConsensusClusterPlus/watchers"
     aria-label="1 user is watching this repository">
    1
  </a>

  </li>

  <li>
        <a href="/login?return_to=%2Frenzhonglu%2FConsensusClusterPlus"
      class="btn btn-sm btn-with-count tooltipped tooltipped-s"
      aria-label="You must be signed in to star a repository" rel="nofollow">
      <svg class="octicon octicon-star v-align-text-bottom" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M14 6l-4.9-.64L7 1 4.9 5.36 0 6l3.6 3.26L2.67 14 7 11.67 11.33 14l-.93-4.74L14 6z"/></svg>
      Star
    </a>

    <a class="social-count js-social-count" href="/renzhonglu/ConsensusClusterPlus/stargazers"
      aria-label="1 user starred this repository">
      1
    </a>

  </li>

  <li>
      <a href="/login?return_to=%2Frenzhonglu%2FConsensusClusterPlus"
        class="btn btn-sm btn-with-count tooltipped tooltipped-s"
        aria-label="You must be signed in to fork a repository" rel="nofollow">
        <svg class="octicon octicon-repo-forked v-align-text-bottom" viewBox="0 0 10 16" version="1.1" width="10" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8 1a1.993 1.993 0 0 0-1 3.72V6L5 8 3 6V4.72A1.993 1.993 0 0 0 2 1a1.993 1.993 0 0 0-1 3.72V6.5l3 3v1.78A1.993 1.993 0 0 0 5 15a1.993 1.993 0 0 0 1-3.72V9.5l3-3V4.72A1.993 1.993 0 0 0 8 1zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3 10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zm3-10c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
        Fork
      </a>

    <a href="/renzhonglu/ConsensusClusterPlus/network/members" class="social-count"
       aria-label="4 users forked this repository">
      4
    </a>
  </li>
</ul>

      <h1 class="public ">
  <svg class="octicon octicon-repo" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M4 9H3V8h1v1zm0-3H3v1h1V6zm0-2H3v1h1V4zm0-2H3v1h1V2zm8-1v12c0 .55-.45 1-1 1H6v2l-1.5-1.5L3 16v-2H1c-.55 0-1-.45-1-1V1c0-.55.45-1 1-1h10c.55 0 1 .45 1 1zm-1 10H1v2h2v-1h3v1h5v-2zm0-10H2v9h9V1z"/></svg>
  <span class="author" itemprop="author"><a class="url fn" rel="author" data-hovercard-type="user" data-hovercard-url="/hovercards?user_id=8721120" data-octo-click="hovercard-link-click" data-octo-dimensions="link_type:self" href="/renzhonglu">renzhonglu</a></span><!--
--><span class="path-divider">/</span><!--
--><strong itemprop="name"><a data-pjax="#js-repo-pjax-container" href="/renzhonglu/ConsensusClusterPlus">ConsensusClusterPlus</a></strong>

</h1>

    </div>
    
<nav class="reponav js-repo-nav js-sidenav-container-pjax container"
     itemscope
     itemtype="http://schema.org/BreadcrumbList"
    aria-label="Repository"
     data-pjax="#js-repo-pjax-container">

  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a class="js-selected-navigation-item selected reponav-item" itemprop="url" data-hotkey="g c" aria-current="page" data-selected-links="repo_source repo_downloads repo_commits repo_releases repo_tags repo_branches repo_packages /renzhonglu/ConsensusClusterPlus" href="/renzhonglu/ConsensusClusterPlus">
      <svg class="octicon octicon-code" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M9.5 3L8 4.5 11.5 8 8 11.5 9.5 13 14 8 9.5 3zm-5 0L0 8l4.5 5L6 11.5 2.5 8 6 4.5 4.5 3z"/></svg>
      <span itemprop="name">Code</span>
      <meta itemprop="position" content="1">
</a>  </span>


  <span itemscope itemtype="http://schema.org/ListItem" itemprop="itemListElement">
    <a data-hotkey="g p" itemprop="url" class="js-selected-navigation-item reponav-item" data-selected-links="repo_pulls checks /renzhonglu/ConsensusClusterPlus/pulls" href="/renzhonglu/ConsensusClusterPlus/pulls">
      <svg class="octicon octicon-git-pull-request" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M11 11.28V5c-.03-.78-.34-1.47-.94-2.06C9.46 2.35 8.78 2.03 8 2H7V0L4 3l3 3V4h1c.27.02.48.11.69.31.21.2.3.42.31.69v6.28A1.993 1.993 0 0 0 10 15a1.993 1.993 0 0 0 1-3.72zm-1 2.92c-.66 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2zM4 3c0-1.11-.89-2-2-2a1.993 1.993 0 0 0-1 3.72v6.56A1.993 1.993 0 0 0 2 15a1.993 1.993 0 0 0 1-3.72V4.72c.59-.34 1-.98 1-1.72zm-.8 10c0 .66-.55 1.2-1.2 1.2-.65 0-1.2-.55-1.2-1.2 0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2zM2 4.2C1.34 4.2.8 3.65.8 3c0-.65.55-1.2 1.2-1.2.65 0 1.2.55 1.2 1.2 0 .65-.55 1.2-1.2 1.2z"/></svg>
      <span itemprop="name">Pull requests</span>
      <span class="Counter">0</span>
      <meta itemprop="position" content="3">
</a>  </span>


    <a data-hotkey="g b" class="js-selected-navigation-item reponav-item" data-selected-links="repo_projects new_repo_project repo_project /renzhonglu/ConsensusClusterPlus/projects" href="/renzhonglu/ConsensusClusterPlus/projects">
      <svg class="octicon octicon-project" viewBox="0 0 15 16" version="1.1" width="15" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M10 12h3V2h-3v10zm-4-2h3V2H6v8zm-4 4h3V2H2v12zm-1 1h13V1H1v14zM14 0H1a1 1 0 0 0-1 1v14a1 1 0 0 0 1 1h13a1 1 0 0 0 1-1V1a1 1 0 0 0-1-1z"/></svg>
      Projects
      <span class="Counter" >0</span>
</a>


    <a class="js-selected-navigation-item reponav-item" data-selected-links="repo_graphs repo_contributors dependency_graph pulse alerts security people /renzhonglu/ConsensusClusterPlus/pulse" href="/renzhonglu/ConsensusClusterPlus/pulse">
      <svg class="octicon octicon-graph" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M16 14v1H0V0h1v14h15zM5 13H3V8h2v5zm4 0H7V3h2v10zm4 0h-2V6h2v7z"/></svg>
      Insights
</a>

</nav>


  </div>
<div class="container new-discussion-timeline experiment-repo-nav  ">
  <div class="repository-content ">

    
    



  
    <a class="d-none js-permalink-shortcut" data-hotkey="y" href="/renzhonglu/ConsensusClusterPlus/blob/08bb6aa1e93619807738d015d8c40f4bf301266b/R/ConsensusClusterPlus.R">Permalink</a>

    <!-- blob contrib key: blob_contributors:v21:a2560aa0de21493c725e319097b16fb4 -->

        <div class="signup-prompt-bg rounded-1">
      <div class="signup-prompt p-4 text-center mb-4 rounded-1">
        <div class="position-relative">
          <!-- '"` --><!-- </textarea></xmp> --></option></form><form action="/site/dismiss_signup_prompt" accept-charset="UTF-8" method="post"><input name="utf8" type="hidden" value="&#x2713;" /><input type="hidden" name="authenticity_token" value="CDhXF4SKosEwCRa1mH19o3QsO8OYInwSjIpwTo6Ly2fAd12NHzx0aTr05lPN0YcN8IaNWmXLdIJ11raDBVqc0g==" />
            <button type="submit" class="position-absolute top-0 right-0 btn-link link-gray" data-ga-click="(Logged out) Sign up prompt, clicked Dismiss, text:dismiss">
              Dismiss
            </button>
</form>          <h3 class="pt-2">Join GitHub today</h3>
          <p class="col-6 mx-auto">GitHub is home to over 31 million developers working together to host and review code, manage projects, and build software together.</p>
          <a class="btn btn-primary" href="/join?source=prompt-blob-show" data-ga-click="(Logged out) Sign up prompt, clicked Sign up, text:sign-up">Sign up</a>
        </div>
      </div>
    </div>


    <div class="file-navigation">
      
<details class="details-reset details-overlay select-menu branch-select-menu float-left">
  <summary class="btn btn-sm select-menu-button css-truncate"
           data-hotkey="w"
           
           title="Switch branches or tags">
    <i>Branch:</i>
    <span class="css-truncate-target">master</span>
  </summary>

  <details-menu class="select-menu-modal position-absolute" style="z-index: 99;" src="/renzhonglu/ConsensusClusterPlus/ref-list/master/R/ConsensusClusterPlus.R?source_action=show&amp;source_controller=blob" preload>
    <include-fragment class="select-menu-loading-overlay anim-pulse">
      <svg height="32" class="octicon octicon-octoface" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M14.7 5.34c.13-.32.55-1.59-.13-3.31 0 0-1.05-.33-3.44 1.3-1-.28-2.07-.32-3.13-.32s-2.13.04-3.13.32c-2.39-1.64-3.44-1.3-3.44-1.3-.68 1.72-.26 2.99-.13 3.31C.49 6.21 0 7.33 0 8.69 0 13.84 3.33 15 7.98 15S16 13.84 16 8.69c0-1.36-.49-2.48-1.3-3.35zM8 14.02c-3.3 0-5.98-.15-5.98-3.35 0-.76.38-1.48 1.02-2.07 1.07-.98 2.9-.46 4.96-.46 2.07 0 3.88-.52 4.96.46.65.59 1.02 1.3 1.02 2.07 0 3.19-2.68 3.35-5.98 3.35zM5.49 9.01c-.66 0-1.2.8-1.2 1.78s.54 1.79 1.2 1.79c.66 0 1.2-.8 1.2-1.79s-.54-1.78-1.2-1.78zm5.02 0c-.66 0-1.2.79-1.2 1.78s.54 1.79 1.2 1.79c.66 0 1.2-.8 1.2-1.79s-.53-1.78-1.2-1.78z"/></svg>
    </include-fragment>
  </details-menu>
</details>

      <div class="BtnGroup float-right">
        <a href="/renzhonglu/ConsensusClusterPlus/find/master"
              class="js-pjax-capture-input btn btn-sm BtnGroup-item"
              data-pjax
              data-hotkey="t">
          Find file
        </a>
        <clipboard-copy for="blob-path" class="btn btn-sm BtnGroup-item">
          Copy path
        </clipboard-copy>
      </div>
      <div id="blob-path" class="breadcrumb">
        <span class="repo-root js-repo-root"><span class="js-path-segment"><a data-pjax="true" href="/renzhonglu/ConsensusClusterPlus"><span>ConsensusClusterPlus</span></a></span></span><span class="separator">/</span><span class="js-path-segment"><a data-pjax="true" href="/renzhonglu/ConsensusClusterPlus/tree/master/R"><span>R</span></a></span><span class="separator">/</span><strong class="final-path">ConsensusClusterPlus.R</strong>
      </div>
    </div>



    
  <div class="commit-tease d-flex flex-column flex-shrink-0">
      <div class="d-flex flex-justify-between ">
        <span class="pr-md-4">
          <img class="avatar" width="20" height="20" alt="" src="https://camo.githubusercontent.com/8786ea7a8091399e16043f8f4ea1eb8365cf89b2/68747470733a2f2f322e67726176617461722e636f6d2f6176617461722f32666366643664636466346539336337613930623135343663666636373661623f643d68747470732533412532462532466769746875622e6769746875626173736574732e636f6d253246696d6167657325324667726176617461727325324667726176617461722d757365722d3432302e706e6726723d6726733d313430" data-canonical-src="https://2.gravatar.com/avatar/2fcfd6dcdf4e93c7a90b1546cff676ab?d=https%3A%2F%2Fgithub.githubassets.com%2Fimages%2Fgravatars%2Fgravatar-user-420.png&amp;r=g&amp;s=140" />
          <span class="user-mention">m.wilkerson</span>
            <a data-pjax="true" title="pFeature less than 1 fix for correlation distance

git-svn-id: https://hedgehog.fhcrc.org/bioconductor/trunk/madman/Rpacks/ConsensusClusterPlus@79900 bc3139a8-67e5-0310-9ffc-ced21a209358" class="message" href="/renzhonglu/ConsensusClusterPlus/commit/77f83da386e16361c8f0919c48fd6e734ca3c79e">pFeature less than 1 fix for correlation distance</a>
        </span>
        <span class="d-inline-block flex-shrink-0 v-align-bottom ">
          <a class="commit-tease-sha pr-2" href="/renzhonglu/ConsensusClusterPlus/commit/77f83da386e16361c8f0919c48fd6e734ca3c79e" data-pjax>
            77f83da
          </a>
          <relative-time datetime="2013-08-28T20:50:34Z">Aug 28, 2013</relative-time>
        </span>
      </div>

    <div class="commit-tease-contributors flex-auto">
      
<details class="details-reset details-overlay details-overlay-dark lh-default text-gray-dark float-left mr-2" id="blob_contributors_box">
  <summary
      class="btn-link"
      aria-haspopup="dialog"
      
      
      >
    
    <span><strong>0</strong> contributors</span>
  </summary>
  <details-dialog class="Box Box--overlay d-flex flex-column anim-fade-in fast " aria-label="Users who have contributed to this file">
    <div class="Box-header">
      <button class="Box-btn-octicon btn-octicon float-right" type="button" aria-label="Close dialog" data-close-dialog>
        <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
      </button>
      <h3 class="Box-title">Users who have contributed to this file</h3>
    </div>
    
        <ul class="list-style-none overflow-auto">
        </ul>

  </details-dialog>
</details>
      
    </div>
  </div>





    <div class="file ">
      
<div class="file-header ">

  <div class="file-info float-left ">
      670 lines (590 sloc)
      <span class="file-info-divider"></span>
    24.5 KB
  </div>

  <div class="file-actions d-flex ">

    <div class="BtnGroup">
      <a id="raw-url" class="btn btn-sm BtnGroup-item" href="/renzhonglu/ConsensusClusterPlus/raw/master/R/ConsensusClusterPlus.R">Raw</a>
        <a class="btn btn-sm js-update-url-with-hash BtnGroup-item" data-hotkey="b" href="/renzhonglu/ConsensusClusterPlus/blame/master/R/ConsensusClusterPlus.R">Blame</a>
      <a rel="nofollow" class="btn btn-sm BtnGroup-item" href="/renzhonglu/ConsensusClusterPlus/commits/master/R/ConsensusClusterPlus.R">History</a>
    </div>


    <div>

          <button type="button" class="btn-octicon disabled tooltipped tooltipped-nw"
            aria-label="You must be signed in to make or propose changes">
            <svg class="octicon octicon-pencil" viewBox="0 0 14 16" version="1.1" width="14" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M0 12v3h3l8-8-3-3-8 8zm3 2H1v-2h1v1h1v1zm10.3-9.3L12 6 9 3l1.3-1.3a.996.996 0 0 1 1.41 0l1.59 1.59c.39.39.39 1.02 0 1.41z"/></svg>
          </button>
          <button type="button" class="btn-octicon btn-octicon-danger disabled tooltipped tooltipped-nw"
            aria-label="You must be signed in to make or propose changes">
            <svg class="octicon octicon-trashcan" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M11 2H9c0-.55-.45-1-1-1H5c-.55 0-1 .45-1 1H2c-.55 0-1 .45-1 1v1c0 .55.45 1 1 1v9c0 .55.45 1 1 1h7c.55 0 1-.45 1-1V5c.55 0 1-.45 1-1V3c0-.55-.45-1-1-1zm-1 12H3V5h1v8h1V5h1v8h1V5h1v8h1V5h1v9zm1-10H2V3h9v1z"/></svg>
          </button>
    </div>
  </div>
</div>

      

  <div itemprop="text" class="blob-wrapper data type-r ">
      
<table class="highlight tab-size js-file-line-container" data-tab-size="8">
      <tr>
        <td id="L1" class="blob-num js-line-number" data-line-number="1"></td>
        <td id="LC1" class="blob-code blob-code-inner js-file-line"><span class="pl-en">ConsensusClusterPlus</span> <span class="pl-k">&lt;-</span> <span class="pl-k">function</span>( <span class="pl-v">d</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L2" class="blob-num js-line-number" data-line-number="2"></td>
        <td id="LC2" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">maxK</span> <span class="pl-k">=</span> <span class="pl-c1">3</span>,</td>
      </tr>
      <tr>
        <td id="L3" class="blob-num js-line-number" data-line-number="3"></td>
        <td id="LC3" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">reps</span><span class="pl-k">=</span><span class="pl-c1">10</span>,</td>
      </tr>
      <tr>
        <td id="L4" class="blob-num js-line-number" data-line-number="4"></td>
        <td id="LC4" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">pItem</span><span class="pl-k">=</span><span class="pl-c1">0.8</span>,</td>
      </tr>
      <tr>
        <td id="L5" class="blob-num js-line-number" data-line-number="5"></td>
        <td id="LC5" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">pFeature</span><span class="pl-k">=</span><span class="pl-c1">1</span>,</td>
      </tr>
      <tr>
        <td id="L6" class="blob-num js-line-number" data-line-number="6"></td>
        <td id="LC6" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">clusterAlg</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>hc<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L7" class="blob-num js-line-number" data-line-number="7"></td>
        <td id="LC7" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>untitled_consensus_cluster<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L8" class="blob-num js-line-number" data-line-number="8"></td>
        <td id="LC8" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">innerLinkage</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>average<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L9" class="blob-num js-line-number" data-line-number="9"></td>
        <td id="LC9" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">finalLinkage</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>average<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L10" class="blob-num js-line-number" data-line-number="10"></td>
        <td id="LC10" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">distance</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>pearson<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L11" class="blob-num js-line-number" data-line-number="11"></td>
        <td id="LC11" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">ml</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L12" class="blob-num js-line-number" data-line-number="12"></td>
        <td id="LC12" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">tmyPal</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L13" class="blob-num js-line-number" data-line-number="13"></td>
        <td id="LC13" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">seed</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L14" class="blob-num js-line-number" data-line-number="14"></td>
        <td id="LC14" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">plot</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L15" class="blob-num js-line-number" data-line-number="15"></td>
        <td id="LC15" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">writeTable</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span>,</td>
      </tr>
      <tr>
        <td id="L16" class="blob-num js-line-number" data-line-number="16"></td>
        <td id="LC16" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">weightsItem</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L17" class="blob-num js-line-number" data-line-number="17"></td>
        <td id="LC17" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">weightsFeature</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L18" class="blob-num js-line-number" data-line-number="18"></td>
        <td id="LC18" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-v">verbose</span><span class="pl-k">=</span><span class="pl-c1">F</span>,</td>
      </tr>
      <tr>
        <td id="L19" class="blob-num js-line-number" data-line-number="19"></td>
        <td id="LC19" class="blob-code blob-code-inner js-file-line">				  <span class="pl-v">corUse</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>everything<span class="pl-pds">&quot;</span></span> ) {</td>
      </tr>
      <tr>
        <td id="L20" class="blob-num js-line-number" data-line-number="20"></td>
        <td id="LC20" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#description: runs consensus subsamples </span></td>
      </tr>
      <tr>
        <td id="L21" class="blob-num js-line-number" data-line-number="21"></td>
        <td id="LC21" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">seed</span>)<span class="pl-k">==</span><span class="pl-c1">TRUE</span>){</td>
      </tr>
      <tr>
        <td id="L22" class="blob-num js-line-number" data-line-number="22"></td>
        <td id="LC22" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">seed</span><span class="pl-k">=</span><span class="pl-v">timeSeed</span> <span class="pl-k">=</span> as.numeric(Sys.time())</td>
      </tr>
      <tr>
        <td id="L23" class="blob-num js-line-number" data-line-number="23"></td>
        <td id="LC23" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L24" class="blob-num js-line-number" data-line-number="24"></td>
        <td id="LC24" class="blob-code blob-code-inner js-file-line">  set.seed(<span class="pl-smi">seed</span>)</td>
      </tr>
      <tr>
        <td id="L25" class="blob-num js-line-number" data-line-number="25"></td>
        <td id="LC25" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L26" class="blob-num js-line-number" data-line-number="26"></td>
        <td id="LC26" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>distance=ifelse( inherits(d,&quot;dist&quot;), attr( d, &quot;method&quot; ), &quot;pearson&quot; )</span></td>
      </tr>
      <tr>
        <td id="L27" class="blob-num js-line-number" data-line-number="27"></td>
        <td id="LC27" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L28" class="blob-num js-line-number" data-line-number="28"></td>
        <td id="LC28" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L29" class="blob-num js-line-number" data-line-number="29"></td>
        <td id="LC29" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">ml</span>)<span class="pl-k">==</span><span class="pl-c1">TRUE</span>){</td>
      </tr>
      <tr>
        <td id="L30" class="blob-num js-line-number" data-line-number="30"></td>
        <td id="LC30" class="blob-code blob-code-inner js-file-line">    </td>
      </tr>
      <tr>
        <td id="L31" class="blob-num js-line-number" data-line-number="31"></td>
        <td id="LC31" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( <span class="pl-k">!</span> class( <span class="pl-smi">d</span> ) <span class="pl-k">%in%</span> c( <span class="pl-s"><span class="pl-pds">&quot;</span>dist<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>matrix<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>ExpressionSet<span class="pl-pds">&quot;</span></span> ) ) {</td>
      </tr>
      <tr>
        <td id="L32" class="blob-num js-line-number" data-line-number="32"></td>
        <td id="LC32" class="blob-code blob-code-inner js-file-line">      stop(<span class="pl-s"><span class="pl-pds">&quot;</span>d must be a matrix, distance object or ExpressionSet (eset object)<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L33" class="blob-num js-line-number" data-line-number="33"></td>
        <td id="LC33" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L34" class="blob-num js-line-number" data-line-number="34"></td>
        <td id="LC34" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L35" class="blob-num js-line-number" data-line-number="35"></td>
        <td id="LC35" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( inherits( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>dist<span class="pl-pds">&quot;</span></span> ) ) {</td>
      </tr>
      <tr>
        <td id="L36" class="blob-num js-line-number" data-line-number="36"></td>
        <td id="LC36" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># if d is a distance matrix, fix a few things so that they don&#39;t cause problems with the analysis</span></td>
      </tr>
      <tr>
        <td id="L37" class="blob-num js-line-number" data-line-number="37"></td>
        <td id="LC37" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#  Note, assumption is that if d is a distance matrix, the user doesn&#39;t want to sample over the row features</span></td>
      </tr>
      <tr>
        <td id="L38" class="blob-num js-line-number" data-line-number="38"></td>
        <td id="LC38" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( is.null( attr( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) ) ) {</td>
      </tr>
      <tr>
        <td id="L39" class="blob-num js-line-number" data-line-number="39"></td>
        <td id="LC39" class="blob-code blob-code-inner js-file-line">        attr( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&lt;-</span> <span class="pl-smi">distance</span> <span class="pl-k">&lt;-</span> <span class="pl-s"><span class="pl-pds">&quot;</span>unknown - user-specified<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L40" class="blob-num js-line-number" data-line-number="40"></td>
        <td id="LC40" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L41" class="blob-num js-line-number" data-line-number="41"></td>
        <td id="LC41" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( is.null( <span class="pl-smi">distance</span> ) <span class="pl-k">||</span> ( <span class="pl-smi">distance</span> <span class="pl-k">!=</span> attr( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) ) ) {</td>
      </tr>
      <tr>
        <td id="L42" class="blob-num js-line-number" data-line-number="42"></td>
        <td id="LC42" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">distance</span> <span class="pl-k">&lt;-</span> attr( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) </td>
      </tr>
      <tr>
        <td id="L43" class="blob-num js-line-number" data-line-number="43"></td>
        <td id="LC43" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L44" class="blob-num js-line-number" data-line-number="44"></td>
        <td id="LC44" class="blob-code blob-code-inner js-file-line">      </td>
      </tr>
      <tr>
        <td id="L45" class="blob-num js-line-number" data-line-number="45"></td>
        <td id="LC45" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( ( <span class="pl-k">!</span> is.null( <span class="pl-smi">pFeature</span> ) ) <span class="pl-k">&amp;&amp;</span> ( <span class="pl-smi">pFeature</span> <span class="pl-k">&lt;</span> <span class="pl-c1">1</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L46" class="blob-num js-line-number" data-line-number="46"></td>
        <td id="LC46" class="blob-code blob-code-inner js-file-line">        message( <span class="pl-s"><span class="pl-pds">&quot;</span>Cannot use the pFeatures parameter when specifying a distance matrix as the data object<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L47" class="blob-num js-line-number" data-line-number="47"></td>
        <td id="LC47" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">pFeature</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">1</span></td>
      </tr>
      <tr>
        <td id="L48" class="blob-num js-line-number" data-line-number="48"></td>
        <td id="LC48" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L49" class="blob-num js-line-number" data-line-number="49"></td>
        <td id="LC49" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( <span class="pl-k">!</span> is.null( <span class="pl-smi">weightsFeature</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L50" class="blob-num js-line-number" data-line-number="50"></td>
        <td id="LC50" class="blob-code blob-code-inner js-file-line">        message( <span class="pl-s"><span class="pl-pds">&quot;</span>Cannot use the weightsFeature parameter when specifying a distance matrix as the data object<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L51" class="blob-num js-line-number" data-line-number="51"></td>
        <td id="LC51" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">weightsFeature</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">NULL</span></td>
      </tr>
      <tr>
        <td id="L52" class="blob-num js-line-number" data-line-number="52"></td>
        <td id="LC52" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L53" class="blob-num js-line-number" data-line-number="53"></td>
        <td id="LC53" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( <span class="pl-smi">clusterAlg</span> <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>km<span class="pl-pds">&quot;</span></span> ) {</td>
      </tr>
      <tr>
        <td id="L54" class="blob-num js-line-number" data-line-number="54"></td>
        <td id="LC54" class="blob-code blob-code-inner js-file-line">        message( <span class="pl-s"><span class="pl-pds">&quot;</span>Note: k-means will cluster the distance matrix you provided.  This is similar to kmdist option when suppling a data matrix<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L55" class="blob-num js-line-number" data-line-number="55"></td>
        <td id="LC55" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span>#d &lt;- as.matrix( d )  #this is now done w/in ccRun</span></td>
      </tr>
      <tr>
        <td id="L56" class="blob-num js-line-number" data-line-number="56"></td>
        <td id="LC56" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L57" class="blob-num js-line-number" data-line-number="57"></td>
        <td id="LC57" class="blob-code blob-code-inner js-file-line">    } <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L58" class="blob-num js-line-number" data-line-number="58"></td>
        <td id="LC58" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( is.null( <span class="pl-smi">distance</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L59" class="blob-num js-line-number" data-line-number="59"></td>
        <td id="LC59" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span># we should never get here, but just in case</span></td>
      </tr>
      <tr>
        <td id="L60" class="blob-num js-line-number" data-line-number="60"></td>
        <td id="LC60" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">distance</span> <span class="pl-k">&lt;-</span> <span class="pl-s"><span class="pl-pds">&quot;</span>pearson<span class="pl-pds">&quot;</span></span></td>
      </tr>
      <tr>
        <td id="L61" class="blob-num js-line-number" data-line-number="61"></td>
        <td id="LC61" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L62" class="blob-num js-line-number" data-line-number="62"></td>
        <td id="LC62" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L63" class="blob-num js-line-number" data-line-number="63"></td>
        <td id="LC63" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L64" class="blob-num js-line-number" data-line-number="64"></td>
        <td id="LC64" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( ( <span class="pl-smi">clusterAlg</span> <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>km<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&amp;&amp;</span> inherits( <span class="pl-smi">distance</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>character<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&amp;&amp;</span> ( <span class="pl-smi">distance</span> <span class="pl-k">!=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>euclidean<span class="pl-pds">&quot;</span></span> ) ) {</td>
      </tr>
      <tr>
        <td id="L65" class="blob-num js-line-number" data-line-number="65"></td>
        <td id="LC65" class="blob-code blob-code-inner js-file-line">      message( <span class="pl-s"><span class="pl-pds">&quot;</span>Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix using k-means use the &#39;kmdist&#39; option, or use a different algorithm such as &#39;hc&#39; or &#39;pam&#39;.  Changing distance to euclidean<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L66" class="blob-num js-line-number" data-line-number="66"></td>
        <td id="LC66" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">distance</span> <span class="pl-k">&lt;-</span> <span class="pl-s"><span class="pl-pds">&#39;</span>euclidean<span class="pl-pds">&#39;</span></span></td>
      </tr>
      <tr>
        <td id="L67" class="blob-num js-line-number" data-line-number="67"></td>
        <td id="LC67" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L68" class="blob-num js-line-number" data-line-number="68"></td>
        <td id="LC68" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L69" class="blob-num js-line-number" data-line-number="69"></td>
        <td id="LC69" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L70" class="blob-num js-line-number" data-line-number="70"></td>
        <td id="LC70" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( inherits( <span class="pl-smi">d</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>ExpressionSet<span class="pl-pds">&quot;</span></span> ) ) {</td>
      </tr>
      <tr>
        <td id="L71" class="blob-num js-line-number" data-line-number="71"></td>
        <td id="LC71" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">d</span> <span class="pl-k">&lt;-</span> exprs(<span class="pl-smi">d</span>)</td>
      </tr>
      <tr>
        <td id="L72" class="blob-num js-line-number" data-line-number="72"></td>
        <td id="LC72" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L73" class="blob-num js-line-number" data-line-number="73"></td>
        <td id="LC73" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L74" class="blob-num js-line-number" data-line-number="74"></td>
        <td id="LC74" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">ml</span> <span class="pl-k">&lt;-</span> ccRun( <span class="pl-v">d</span><span class="pl-k">=</span><span class="pl-smi">d</span>,</td>
      </tr>
      <tr>
        <td id="L75" class="blob-num js-line-number" data-line-number="75"></td>
        <td id="LC75" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">maxK</span><span class="pl-k">=</span><span class="pl-smi">maxK</span>,</td>
      </tr>
      <tr>
        <td id="L76" class="blob-num js-line-number" data-line-number="76"></td>
        <td id="LC76" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">repCount</span><span class="pl-k">=</span><span class="pl-smi">reps</span>,</td>
      </tr>
      <tr>
        <td id="L77" class="blob-num js-line-number" data-line-number="77"></td>
        <td id="LC77" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">diss</span><span class="pl-k">=</span>inherits(<span class="pl-smi">d</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>dist<span class="pl-pds">&quot;</span></span>),</td>
      </tr>
      <tr>
        <td id="L78" class="blob-num js-line-number" data-line-number="78"></td>
        <td id="LC78" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">pItem</span><span class="pl-k">=</span><span class="pl-smi">pItem</span>,</td>
      </tr>
      <tr>
        <td id="L79" class="blob-num js-line-number" data-line-number="79"></td>
        <td id="LC79" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">pFeature</span><span class="pl-k">=</span><span class="pl-smi">pFeature</span>,</td>
      </tr>
      <tr>
        <td id="L80" class="blob-num js-line-number" data-line-number="80"></td>
        <td id="LC80" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">innerLinkage</span><span class="pl-k">=</span><span class="pl-smi">innerLinkage</span>,</td>
      </tr>
      <tr>
        <td id="L81" class="blob-num js-line-number" data-line-number="81"></td>
        <td id="LC81" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">clusterAlg</span><span class="pl-k">=</span><span class="pl-smi">clusterAlg</span>,</td>
      </tr>
      <tr>
        <td id="L82" class="blob-num js-line-number" data-line-number="82"></td>
        <td id="LC82" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">weightsFeature</span><span class="pl-k">=</span><span class="pl-smi">weightsFeature</span>,</td>
      </tr>
      <tr>
        <td id="L83" class="blob-num js-line-number" data-line-number="83"></td>
        <td id="LC83" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">weightsItem</span><span class="pl-k">=</span><span class="pl-smi">weightsItem</span>,</td>
      </tr>
      <tr>
        <td id="L84" class="blob-num js-line-number" data-line-number="84"></td>
        <td id="LC84" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">distance</span><span class="pl-k">=</span><span class="pl-smi">distance</span>,</td>
      </tr>
      <tr>
        <td id="L85" class="blob-num js-line-number" data-line-number="85"></td>
        <td id="LC85" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">verbose</span><span class="pl-k">=</span><span class="pl-smi">verbose</span>,</td>
      </tr>
      <tr>
        <td id="L86" class="blob-num js-line-number" data-line-number="86"></td>
        <td id="LC86" class="blob-code blob-code-inner js-file-line">		 <span class="pl-v">corUse</span><span class="pl-k">=</span><span class="pl-smi">corUse</span>)</td>
      </tr>
      <tr>
        <td id="L87" class="blob-num js-line-number" data-line-number="87"></td>
        <td id="LC87" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L88" class="blob-num js-line-number" data-line-number="88"></td>
        <td id="LC88" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">res</span><span class="pl-k">=</span><span class="pl-k">list</span>();</td>
      </tr>
      <tr>
        <td id="L89" class="blob-num js-line-number" data-line-number="89"></td>
        <td id="LC89" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L90" class="blob-num js-line-number" data-line-number="90"></td>
        <td id="LC90" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#make results directory</span></td>
      </tr>
      <tr>
        <td id="L91" class="blob-num js-line-number" data-line-number="91"></td>
        <td id="LC91" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>((is.null(<span class="pl-smi">plot</span>)<span class="pl-k">==</span><span class="pl-c1">FALSE</span> <span class="pl-k">|</span> <span class="pl-smi">writeTable</span>) <span class="pl-k">&amp;</span> <span class="pl-k">!</span>file.exists(paste(<span class="pl-smi">title</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))){</td>
      </tr>
      <tr>
        <td id="L92" class="blob-num js-line-number" data-line-number="92"></td>
        <td id="LC92" class="blob-code blob-code-inner js-file-line">    dir.create(paste(<span class="pl-smi">title</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L93" class="blob-num js-line-number" data-line-number="93"></td>
        <td id="LC93" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L94" class="blob-num js-line-number" data-line-number="94"></td>
        <td id="LC94" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L95" class="blob-num js-line-number" data-line-number="95"></td>
        <td id="LC95" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#write log file</span></td>
      </tr>
      <tr>
        <td id="L96" class="blob-num js-line-number" data-line-number="96"></td>
        <td id="LC96" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">log</span> <span class="pl-k">&lt;-</span> <span class="pl-k">matrix</span>( <span class="pl-v">ncol</span><span class="pl-k">=</span><span class="pl-c1">2</span>,</td>
      </tr>
      <tr>
        <td id="L97" class="blob-num js-line-number" data-line-number="97"></td>
        <td id="LC97" class="blob-code blob-code-inner js-file-line">                 <span class="pl-v">byrow</span><span class="pl-k">=</span><span class="pl-c1">T</span>,</td>
      </tr>
      <tr>
        <td id="L98" class="blob-num js-line-number" data-line-number="98"></td>
        <td id="LC98" class="blob-code blob-code-inner js-file-line">                 c(<span class="pl-s"><span class="pl-pds">&quot;</span>title<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">title</span>,</td>
      </tr>
      <tr>
        <td id="L99" class="blob-num js-line-number" data-line-number="99"></td>
        <td id="LC99" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>maxK<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">maxK</span>,</td>
      </tr>
      <tr>
        <td id="L100" class="blob-num js-line-number" data-line-number="100"></td>
        <td id="LC100" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>input matrix rows<span class="pl-pds">&quot;</span></span>,ifelse ( inherits( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>matrix<span class="pl-pds">&quot;</span></span> ), nrow(<span class="pl-smi">d</span>), <span class="pl-s"><span class="pl-pds">&quot;</span>dist-mat<span class="pl-pds">&quot;</span></span> ), </td>
      </tr>
      <tr>
        <td id="L101" class="blob-num js-line-number" data-line-number="101"></td>
        <td id="LC101" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>input matrix columns<span class="pl-pds">&quot;</span></span>,ifelse ( inherits( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>matrix<span class="pl-pds">&quot;</span></span> ), ncol(<span class="pl-smi">d</span>), ncol( as.matrix(<span class="pl-smi">d</span>) ) ), </td>
      </tr>
      <tr>
        <td id="L102" class="blob-num js-line-number" data-line-number="102"></td>
        <td id="LC102" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>number of bootstraps<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">reps</span>,</td>
      </tr>
      <tr>
        <td id="L103" class="blob-num js-line-number" data-line-number="103"></td>
        <td id="LC103" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>item subsampling proportion<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">pItem</span>,</td>
      </tr>
      <tr>
        <td id="L104" class="blob-num js-line-number" data-line-number="104"></td>
        <td id="LC104" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>feature subsampling proportion<span class="pl-pds">&quot;</span></span>,ifelse( is.null(<span class="pl-smi">pFeature</span>), <span class="pl-c1">1</span>, <span class="pl-smi">pFeature</span> ),</td>
      </tr>
      <tr>
        <td id="L105" class="blob-num js-line-number" data-line-number="105"></td>
        <td id="LC105" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>cluster algorithm<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">clusterAlg</span>,</td>
      </tr>
      <tr>
        <td id="L106" class="blob-num js-line-number" data-line-number="106"></td>
        <td id="LC106" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>inner linkage type<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">innerLinkage</span>,</td>
      </tr>
      <tr>
        <td id="L107" class="blob-num js-line-number" data-line-number="107"></td>
        <td id="LC107" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>final linkage type<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">finalLinkage</span>,</td>
      </tr>
      <tr>
        <td id="L108" class="blob-num js-line-number" data-line-number="108"></td>
        <td id="LC108" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>correlation method<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">distance</span>,</td>
      </tr>
      <tr>
        <td id="L109" class="blob-num js-line-number" data-line-number="109"></td>
        <td id="LC109" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>plot<span class="pl-pds">&quot;</span></span>,<span class="pl-k">if</span>(is.null(<span class="pl-smi">plot</span>)) <span class="pl-c1">NA</span> <span class="pl-k">else</span> <span class="pl-smi">plot</span>,</td>
      </tr>
      <tr>
        <td id="L110" class="blob-num js-line-number" data-line-number="110"></td>
        <td id="LC110" class="blob-code blob-code-inner js-file-line">                   <span class="pl-s"><span class="pl-pds">&quot;</span>seed<span class="pl-pds">&quot;</span></span>,<span class="pl-k">if</span>(is.null(<span class="pl-smi">seed</span>)) <span class="pl-c1">NA</span> <span class="pl-k">else</span> <span class="pl-smi">seed</span>))</td>
      </tr>
      <tr>
        <td id="L111" class="blob-num js-line-number" data-line-number="111"></td>
        <td id="LC111" class="blob-code blob-code-inner js-file-line">  colnames(<span class="pl-smi">log</span>) <span class="pl-k">=</span> c(<span class="pl-s"><span class="pl-pds">&quot;</span>argument<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>value<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L112" class="blob-num js-line-number" data-line-number="112"></td>
        <td id="LC112" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(<span class="pl-smi">writeTable</span>){</td>
      </tr>
      <tr>
        <td id="L113" class="blob-num js-line-number" data-line-number="113"></td>
        <td id="LC113" class="blob-code blob-code-inner js-file-line">    write.csv(<span class="pl-v">file</span><span class="pl-k">=</span>paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.log.csv<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>), <span class="pl-smi">log</span>,<span class="pl-v">row.names</span><span class="pl-k">=</span><span class="pl-c1">F</span>)</td>
      </tr>
      <tr>
        <td id="L114" class="blob-num js-line-number" data-line-number="114"></td>
        <td id="LC114" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L115" class="blob-num js-line-number" data-line-number="115"></td>
        <td id="LC115" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">plot</span>)){</td>
      </tr>
      <tr>
        <td id="L116" class="blob-num js-line-number" data-line-number="116"></td>
        <td id="LC116" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>#nothing</span></td>
      </tr>
      <tr>
        <td id="L117" class="blob-num js-line-number" data-line-number="117"></td>
        <td id="LC117" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>pngBMP<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L118" class="blob-num js-line-number" data-line-number="118"></td>
        <td id="LC118" class="blob-code blob-code-inner js-file-line">    bitmap(paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>consensus%03d.png<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L119" class="blob-num js-line-number" data-line-number="119"></td>
        <td id="LC119" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>png<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L120" class="blob-num js-line-number" data-line-number="120"></td>
        <td id="LC120" class="blob-code blob-code-inner js-file-line">    png(paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>consensus%03d.png<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L121" class="blob-num js-line-number" data-line-number="121"></td>
        <td id="LC121" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L122" class="blob-num js-line-number" data-line-number="122"></td>
        <td id="LC122" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span> (<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>pdf<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L123" class="blob-num js-line-number" data-line-number="123"></td>
        <td id="LC123" class="blob-code blob-code-inner js-file-line">    pdf(<span class="pl-v">onefile</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>, paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>consensus.pdf<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L124" class="blob-num js-line-number" data-line-number="124"></td>
        <td id="LC124" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span> (<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>ps<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L125" class="blob-num js-line-number" data-line-number="125"></td>
        <td id="LC125" class="blob-code blob-code-inner js-file-line">    postscript(<span class="pl-v">onefile</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>, paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>consensus.ps<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L126" class="blob-num js-line-number" data-line-number="126"></td>
        <td id="LC126" class="blob-code blob-code-inner js-file-line">  }	</td>
      </tr>
      <tr>
        <td id="L127" class="blob-num js-line-number" data-line-number="127"></td>
        <td id="LC127" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L128" class="blob-num js-line-number" data-line-number="128"></td>
        <td id="LC128" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">colorList</span><span class="pl-k">=</span><span class="pl-k">list</span>()</td>
      </tr>
      <tr>
        <td id="L129" class="blob-num js-line-number" data-line-number="129"></td>
        <td id="LC129" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">colorM</span> <span class="pl-k">=</span> rbind() <span class="pl-c"><span class="pl-c">#</span>matrix of colors.</span></td>
      </tr>
      <tr>
        <td id="L130" class="blob-num js-line-number" data-line-number="130"></td>
        <td id="LC130" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L131" class="blob-num js-line-number" data-line-number="131"></td>
        <td id="LC131" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>18 colors for marking different clusters</span></td>
      </tr>
      <tr>
        <td id="L132" class="blob-num js-line-number" data-line-number="132"></td>
        <td id="LC132" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">thisPal</span> <span class="pl-k">&lt;-</span> c(<span class="pl-s"><span class="pl-pds">&quot;</span>#A6CEE3<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#1F78B4<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#B2DF8A<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#33A02C<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#FB9A99<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#E31A1C<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#FDBF6F<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#FF7F00<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#CAB2D6<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#6A3D9A<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#FFFF99<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>#B15928<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L133" class="blob-num js-line-number" data-line-number="133"></td>
        <td id="LC133" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#bd18ea<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>magenta</span></td>
      </tr>
      <tr>
        <td id="L134" class="blob-num js-line-number" data-line-number="134"></td>
        <td id="LC134" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#2ef4ca<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>aqua</span></td>
      </tr>
      <tr>
        <td id="L135" class="blob-num js-line-number" data-line-number="135"></td>
        <td id="LC135" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#f4cced<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>pink,</span></td>
      </tr>
      <tr>
        <td id="L136" class="blob-num js-line-number" data-line-number="136"></td>
        <td id="LC136" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#f4cc03<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>lightorange</span></td>
      </tr>
      <tr>
        <td id="L137" class="blob-num js-line-number" data-line-number="137"></td>
        <td id="LC137" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#05188a<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>navy,</span></td>
      </tr>
      <tr>
        <td id="L138" class="blob-num js-line-number" data-line-number="138"></td>
        <td id="LC138" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#e5a25a<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>light brown</span></td>
      </tr>
      <tr>
        <td id="L139" class="blob-num js-line-number" data-line-number="139"></td>
        <td id="LC139" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#06f106<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>bright green</span></td>
      </tr>
      <tr>
        <td id="L140" class="blob-num js-line-number" data-line-number="140"></td>
        <td id="LC140" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#85848f<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>med gray</span></td>
      </tr>
      <tr>
        <td id="L141" class="blob-num js-line-number" data-line-number="141"></td>
        <td id="LC141" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#000000<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>black</span></td>
      </tr>
      <tr>
        <td id="L142" class="blob-num js-line-number" data-line-number="142"></td>
        <td id="LC142" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#076f25<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>dark green</span></td>
      </tr>
      <tr>
        <td id="L143" class="blob-num js-line-number" data-line-number="143"></td>
        <td id="LC143" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#93cd7f<span class="pl-pds">&quot;</span></span>,<span class="pl-c"><span class="pl-c">#</span>lime green</span></td>
      </tr>
      <tr>
        <td id="L144" class="blob-num js-line-number" data-line-number="144"></td>
        <td id="LC144" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#4d0776<span class="pl-pds">&quot;</span></span>, <span class="pl-c"><span class="pl-c">#</span>dark purple</span></td>
      </tr>
      <tr>
        <td id="L145" class="blob-num js-line-number" data-line-number="145"></td>
        <td id="LC145" class="blob-code blob-code-inner js-file-line">               <span class="pl-s"><span class="pl-pds">&quot;</span>#ffffff<span class="pl-pds">&quot;</span></span> <span class="pl-c"><span class="pl-c">#</span>white</span></td>
      </tr>
      <tr>
        <td id="L146" class="blob-num js-line-number" data-line-number="146"></td>
        <td id="LC146" class="blob-code blob-code-inner js-file-line">               )</td>
      </tr>
      <tr>
        <td id="L147" class="blob-num js-line-number" data-line-number="147"></td>
        <td id="LC147" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L148" class="blob-num js-line-number" data-line-number="148"></td>
        <td id="LC148" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#plot scale</span></td>
      </tr>
      <tr>
        <td id="L149" class="blob-num js-line-number" data-line-number="149"></td>
        <td id="LC149" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">colBreaks</span><span class="pl-k">=</span><span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L150" class="blob-num js-line-number" data-line-number="150"></td>
        <td id="LC150" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">tmyPal</span>)<span class="pl-k">==</span><span class="pl-c1">TRUE</span>){</td>
      </tr>
      <tr>
        <td id="L151" class="blob-num js-line-number" data-line-number="151"></td>
        <td id="LC151" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">colBreaks</span><span class="pl-k">=</span><span class="pl-c1">10</span></td>
      </tr>
      <tr>
        <td id="L152" class="blob-num js-line-number" data-line-number="152"></td>
        <td id="LC152" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">tmyPal</span> <span class="pl-k">=</span> myPal(<span class="pl-smi">colBreaks</span>)</td>
      </tr>
      <tr>
        <td id="L153" class="blob-num js-line-number" data-line-number="153"></td>
        <td id="LC153" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L154" class="blob-num js-line-number" data-line-number="154"></td>
        <td id="LC154" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">colBreaks</span><span class="pl-k">=</span>length(<span class="pl-smi">tmyPal</span>)</td>
      </tr>
      <tr>
        <td id="L155" class="blob-num js-line-number" data-line-number="155"></td>
        <td id="LC155" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L156" class="blob-num js-line-number" data-line-number="156"></td>
        <td id="LC156" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">sc</span> <span class="pl-k">=</span> cbind(seq(<span class="pl-c1">0</span>,<span class="pl-c1">1</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span><span class="pl-k">/</span>( <span class="pl-smi">colBreaks</span>) )); rownames(<span class="pl-smi">sc</span>) <span class="pl-k">=</span> <span class="pl-smi">sc</span>[,<span class="pl-c1">1</span>]</td>
      </tr>
      <tr>
        <td id="L157" class="blob-num js-line-number" data-line-number="157"></td>
        <td id="LC157" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">sc</span> <span class="pl-k">=</span> cbind(<span class="pl-smi">sc</span>,<span class="pl-smi">sc</span>)</td>
      </tr>
      <tr>
        <td id="L158" class="blob-num js-line-number" data-line-number="158"></td>
        <td id="LC158" class="blob-code blob-code-inner js-file-line">  heatmap(<span class="pl-smi">sc</span>, <span class="pl-v">Colv</span><span class="pl-k">=</span><span class="pl-c1">NA</span>, <span class="pl-v">Rowv</span><span class="pl-k">=</span><span class="pl-c1">NA</span>, <span class="pl-v">symm</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span>, <span class="pl-v">scale</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>none<span class="pl-pds">&#39;</span></span>, <span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">tmyPal</span>, <span class="pl-v">na.rm</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>,<span class="pl-v">labRow</span><span class="pl-k">=</span>rownames(<span class="pl-smi">sc</span>),<span class="pl-v">labCol</span><span class="pl-k">=</span><span class="pl-c1">F</span>,<span class="pl-v">main</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>consensus matrix legend<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L159" class="blob-num js-line-number" data-line-number="159"></td>
        <td id="LC159" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L160" class="blob-num js-line-number" data-line-number="160"></td>
        <td id="LC160" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> (<span class="pl-smi">tk</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span><span class="pl-smi">maxK</span>){</td>
      </tr>
      <tr>
        <td id="L161" class="blob-num js-line-number" data-line-number="161"></td>
        <td id="LC161" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(<span class="pl-smi">verbose</span>){</td>
      </tr>
      <tr>
        <td id="L162" class="blob-num js-line-number" data-line-number="162"></td>
        <td id="LC162" class="blob-code blob-code-inner js-file-line">      message(paste(<span class="pl-s"><span class="pl-pds">&quot;</span>consensus <span class="pl-pds">&quot;</span></span>,<span class="pl-smi">tk</span>))</td>
      </tr>
      <tr>
        <td id="L163" class="blob-num js-line-number" data-line-number="163"></td>
        <td id="LC163" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L164" class="blob-num js-line-number" data-line-number="164"></td>
        <td id="LC164" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">fm</span> <span class="pl-k">=</span> <span class="pl-smi">ml</span>[[<span class="pl-smi">tk</span>]]</td>
      </tr>
      <tr>
        <td id="L165" class="blob-num js-line-number" data-line-number="165"></td>
        <td id="LC165" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">hc</span><span class="pl-k">=</span>hclust( as.dist( <span class="pl-c1">1</span> <span class="pl-k">-</span> <span class="pl-smi">fm</span> ), <span class="pl-v">method</span><span class="pl-k">=</span><span class="pl-smi">finalLinkage</span>);</td>
      </tr>
      <tr>
        <td id="L166" class="blob-num js-line-number" data-line-number="166"></td>
        <td id="LC166" class="blob-code blob-code-inner js-file-line">    message(<span class="pl-s"><span class="pl-pds">&quot;</span>clustered<span class="pl-pds">&quot;</span></span>)	</td>
      </tr>
      <tr>
        <td id="L167" class="blob-num js-line-number" data-line-number="167"></td>
        <td id="LC167" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">ct</span> <span class="pl-k">=</span> cutree(<span class="pl-smi">hc</span>,<span class="pl-smi">tk</span>)</td>
      </tr>
      <tr>
        <td id="L168" class="blob-num js-line-number" data-line-number="168"></td>
        <td id="LC168" class="blob-code blob-code-inner js-file-line">    names(<span class="pl-smi">ct</span>) <span class="pl-k">=</span> colnames(<span class="pl-smi">d</span>)</td>
      </tr>
      <tr>
        <td id="L169" class="blob-num js-line-number" data-line-number="169"></td>
        <td id="LC169" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(class(<span class="pl-smi">d</span>)<span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>dist<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L170" class="blob-num js-line-number" data-line-number="170"></td>
        <td id="LC170" class="blob-code blob-code-inner js-file-line">      names(<span class="pl-smi">ct</span>) <span class="pl-k">=</span> colnames(as.matrix(<span class="pl-smi">d</span>))</td>
      </tr>
      <tr>
        <td id="L171" class="blob-num js-line-number" data-line-number="171"></td>
        <td id="LC171" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L172" class="blob-num js-line-number" data-line-number="172"></td>
        <td id="LC172" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">c</span> <span class="pl-k">=</span> <span class="pl-smi">fm</span></td>
      </tr>
      <tr>
        <td id="L173" class="blob-num js-line-number" data-line-number="173"></td>
        <td id="LC173" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L174" class="blob-num js-line-number" data-line-number="174"></td>
        <td id="LC174" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">colorList</span> <span class="pl-k">=</span> setClusterColors(<span class="pl-smi">res</span>[[<span class="pl-smi">tk</span><span class="pl-k">-</span><span class="pl-c1">1</span>]][[<span class="pl-c1">3</span>]],<span class="pl-smi">ct</span>,<span class="pl-smi">thisPal</span>,<span class="pl-smi">colorList</span>)</td>
      </tr>
      <tr>
        <td id="L175" class="blob-num js-line-number" data-line-number="175"></td>
        <td id="LC175" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">pc</span> <span class="pl-k">=</span> <span class="pl-smi">c</span></td>
      </tr>
      <tr>
        <td id="L176" class="blob-num js-line-number" data-line-number="176"></td>
        <td id="LC176" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">pc</span><span class="pl-k">=</span><span class="pl-smi">pc</span>[<span class="pl-smi">hc</span><span class="pl-k">$</span><span class="pl-smi">order</span>,] <span class="pl-c"><span class="pl-c">#</span>pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.</span></td>
      </tr>
      <tr>
        <td id="L177" class="blob-num js-line-number" data-line-number="177"></td>
        <td id="LC177" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L178" class="blob-num js-line-number" data-line-number="178"></td>
        <td id="LC178" class="blob-code blob-code-inner js-file-line">    </td>
      </tr>
      <tr>
        <td id="L179" class="blob-num js-line-number" data-line-number="179"></td>
        <td id="LC179" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L180" class="blob-num js-line-number" data-line-number="180"></td>
        <td id="LC180" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(<span class="pl-k">!</span>is.null(<span class="pl-smi">plot</span>) <span class="pl-k">&amp;&amp;</span> <span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>pngBMP<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L181" class="blob-num js-line-number" data-line-number="181"></td>
        <td id="LC181" class="blob-code blob-code-inner js-file-line">         <span class="pl-v">pc</span> <span class="pl-k">=</span> <span class="pl-smi">pc</span>[,<span class="pl-smi">hc</span><span class="pl-k">$</span><span class="pl-smi">order</span> ] <span class="pl-c"><span class="pl-c">#</span>mod for no tree</span></td>
      </tr>
      <tr>
        <td id="L182" class="blob-num js-line-number" data-line-number="182"></td>
        <td id="LC182" class="blob-code blob-code-inner js-file-line">        <span class="pl-v">pc</span> <span class="pl-k">=</span> rbind(<span class="pl-smi">pc</span>,<span class="pl-c1">0</span>)</td>
      </tr>
      <tr>
        <td id="L183" class="blob-num js-line-number" data-line-number="183"></td>
        <td id="LC183" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">#</span>no dendrogram if pngBMP</span></td>
      </tr>
      <tr>
        <td id="L184" class="blob-num js-line-number" data-line-number="184"></td>
        <td id="LC184" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">oc</span> <span class="pl-k">=</span> <span class="pl-smi">colorList</span>[[<span class="pl-c1">1</span>]][<span class="pl-smi">hc</span><span class="pl-k">$</span><span class="pl-smi">order</span>] <span class="pl-c"><span class="pl-c">#</span>mod for no tree</span></td>
      </tr>
      <tr>
        <td id="L185" class="blob-num js-line-number" data-line-number="185"></td>
        <td id="LC185" class="blob-code blob-code-inner js-file-line">	heatmap(<span class="pl-smi">pc</span>, <span class="pl-v">Colv</span> <span class="pl-k">=</span> <span class="pl-c1">NA</span>, <span class="pl-v">Rowv</span> <span class="pl-k">=</span> <span class="pl-c1">NA</span>, <span class="pl-v">symm</span> <span class="pl-k">=</span> <span class="pl-c1">FALSE</span>, <span class="pl-v">scale</span> <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>none<span class="pl-pds">&quot;</span></span>, <span class="pl-v">col</span> <span class="pl-k">=</span> <span class="pl-smi">tmyPal</span>, <span class="pl-v">na.rm</span> <span class="pl-k">=</span> <span class="pl-c1">TRUE</span>, <span class="pl-v">labRow</span> <span class="pl-k">=</span> <span class="pl-c1">F</span>, <span class="pl-v">labCol</span> <span class="pl-k">=</span> <span class="pl-c1">F</span>, <span class="pl-v">mar</span> <span class="pl-k">=</span> c(<span class="pl-c1">5</span>, <span class="pl-c1">5</span>), <span class="pl-v">main</span> <span class="pl-k">=</span> paste(<span class="pl-s"><span class="pl-pds">&quot;</span>consensus matrix k=<span class="pl-pds">&quot;</span></span>, </td>
      </tr>
      <tr>
        <td id="L186" class="blob-num js-line-number" data-line-number="186"></td>
        <td id="LC186" class="blob-code blob-code-inner js-file-line">                <span class="pl-smi">tk</span>, <span class="pl-v">sep</span> <span class="pl-k">=</span> <span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>), <span class="pl-v">ColSideCol</span> <span class="pl-k">=</span> <span class="pl-smi">oc</span>)</td>
      </tr>
      <tr>
        <td id="L187" class="blob-num js-line-number" data-line-number="187"></td>
        <td id="LC187" class="blob-code blob-code-inner js-file-line">    }<span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L188" class="blob-num js-line-number" data-line-number="188"></td>
        <td id="LC188" class="blob-code blob-code-inner js-file-line">        <span class="pl-v">pc</span> <span class="pl-k">=</span> rbind(<span class="pl-smi">pc</span>,<span class="pl-c1">0</span>)</td>
      </tr>
      <tr>
        <td id="L189" class="blob-num js-line-number" data-line-number="189"></td>
        <td id="LC189" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>former with tree:</span></td>
      </tr>
      <tr>
        <td id="L190" class="blob-num js-line-number" data-line-number="190"></td>
        <td id="LC190" class="blob-code blob-code-inner js-file-line">    	heatmap(<span class="pl-smi">pc</span>, <span class="pl-v">Colv</span><span class="pl-k">=</span>as.dendrogram(<span class="pl-smi">hc</span>), <span class="pl-v">Rowv</span><span class="pl-k">=</span><span class="pl-c1">NA</span>, <span class="pl-v">symm</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span>, <span class="pl-v">scale</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>none<span class="pl-pds">&#39;</span></span>, <span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">tmyPal</span>, <span class="pl-v">na.rm</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>,<span class="pl-v">labRow</span><span class="pl-k">=</span><span class="pl-c1">F</span>,<span class="pl-v">labCol</span><span class="pl-k">=</span><span class="pl-c1">F</span>,<span class="pl-v">mar</span><span class="pl-k">=</span>c(<span class="pl-c1">5</span>,<span class="pl-c1">5</span>),<span class="pl-v">main</span><span class="pl-k">=</span>paste(<span class="pl-s"><span class="pl-pds">&quot;</span>consensus matrix k=<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">tk</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>) , <span class="pl-v">ColSideCol</span><span class="pl-k">=</span><span class="pl-smi">colorList</span>[[<span class="pl-c1">1</span>]])</td>
      </tr>
      <tr>
        <td id="L191" class="blob-num js-line-number" data-line-number="191"></td>
        <td id="LC191" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L192" class="blob-num js-line-number" data-line-number="192"></td>
        <td id="LC192" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L193" class="blob-num js-line-number" data-line-number="193"></td>
        <td id="LC193" class="blob-code blob-code-inner js-file-line">    legend(<span class="pl-s"><span class="pl-pds">&quot;</span>topright<span class="pl-pds">&quot;</span></span>,<span class="pl-v">legend</span><span class="pl-k">=</span>unique(<span class="pl-smi">ct</span>),<span class="pl-v">fill</span><span class="pl-k">=</span>unique(<span class="pl-smi">colorList</span>[[<span class="pl-c1">1</span>]]),<span class="pl-v">horiz</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span> )</td>
      </tr>
      <tr>
        <td id="L194" class="blob-num js-line-number" data-line-number="194"></td>
        <td id="LC194" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L195" class="blob-num js-line-number" data-line-number="195"></td>
        <td id="LC195" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">res</span>[[<span class="pl-smi">tk</span>]] <span class="pl-k">=</span> <span class="pl-k">list</span>(<span class="pl-v">consensusMatrix</span><span class="pl-k">=</span><span class="pl-smi">c</span>,<span class="pl-v">consensusTree</span><span class="pl-k">=</span><span class="pl-smi">hc</span>,<span class="pl-v">consensusClass</span><span class="pl-k">=</span><span class="pl-smi">ct</span>,<span class="pl-v">ml</span><span class="pl-k">=</span><span class="pl-smi">ml</span>[[<span class="pl-smi">tk</span>]],<span class="pl-v">clrs</span><span class="pl-k">=</span><span class="pl-smi">colorList</span>)</td>
      </tr>
      <tr>
        <td id="L196" class="blob-num js-line-number" data-line-number="196"></td>
        <td id="LC196" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">colorM</span> <span class="pl-k">=</span> rbind(<span class="pl-smi">colorM</span>,<span class="pl-smi">colorList</span>[[<span class="pl-c1">1</span>]]) </td>
      </tr>
      <tr>
        <td id="L197" class="blob-num js-line-number" data-line-number="197"></td>
        <td id="LC197" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L198" class="blob-num js-line-number" data-line-number="198"></td>
        <td id="LC198" class="blob-code blob-code-inner js-file-line">  CDF(<span class="pl-smi">ml</span>)</td>
      </tr>
      <tr>
        <td id="L199" class="blob-num js-line-number" data-line-number="199"></td>
        <td id="LC199" class="blob-code blob-code-inner js-file-line">  clusterTrackingPlot(<span class="pl-smi">colorM</span>[,<span class="pl-smi">res</span>[[length(<span class="pl-smi">res</span>)]]<span class="pl-k">$</span><span class="pl-smi">consensusTree</span><span class="pl-k">$</span><span class="pl-smi">order</span>])</td>
      </tr>
      <tr>
        <td id="L200" class="blob-num js-line-number" data-line-number="200"></td>
        <td id="LC200" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">plot</span>)<span class="pl-k">==</span><span class="pl-c1">FALSE</span>){</td>
      </tr>
      <tr>
        <td id="L201" class="blob-num js-line-number" data-line-number="201"></td>
        <td id="LC201" class="blob-code blob-code-inner js-file-line">    dev.off();</td>
      </tr>
      <tr>
        <td id="L202" class="blob-num js-line-number" data-line-number="202"></td>
        <td id="LC202" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L203" class="blob-num js-line-number" data-line-number="203"></td>
        <td id="LC203" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">res</span>[[<span class="pl-c1">1</span>]] <span class="pl-k">=</span> <span class="pl-smi">colorM</span></td>
      </tr>
      <tr>
        <td id="L204" class="blob-num js-line-number" data-line-number="204"></td>
        <td id="LC204" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(<span class="pl-smi">writeTable</span>){</td>
      </tr>
      <tr>
        <td id="L205" class="blob-num js-line-number" data-line-number="205"></td>
        <td id="LC205" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span>(<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span>length(<span class="pl-smi">res</span>)){</td>
      </tr>
      <tr>
        <td id="L206" class="blob-num js-line-number" data-line-number="206"></td>
        <td id="LC206" class="blob-code blob-code-inner js-file-line">      write.csv(<span class="pl-v">file</span><span class="pl-k">=</span>paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.k=<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">i</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.consensusMatrix.csv<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>), <span class="pl-smi">res</span>[[<span class="pl-smi">i</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusMatrix</span>)</td>
      </tr>
      <tr>
        <td id="L207" class="blob-num js-line-number" data-line-number="207"></td>
        <td id="LC207" class="blob-code blob-code-inner js-file-line">      write.table(<span class="pl-v">file</span><span class="pl-k">=</span>paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.k=<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">i</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.consensusClass.csv<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>), <span class="pl-smi">res</span>[[<span class="pl-smi">i</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusClass</span>,<span class="pl-v">col.names</span> <span class="pl-k">=</span> <span class="pl-c1">F</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>,<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L208" class="blob-num js-line-number" data-line-number="208"></td>
        <td id="LC208" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L209" class="blob-num js-line-number" data-line-number="209"></td>
        <td id="LC209" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L210" class="blob-num js-line-number" data-line-number="210"></td>
        <td id="LC210" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">return</span>(<span class="pl-smi">res</span>)</td>
      </tr>
      <tr>
        <td id="L211" class="blob-num js-line-number" data-line-number="211"></td>
        <td id="LC211" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L212" class="blob-num js-line-number" data-line-number="212"></td>
        <td id="LC212" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L213" class="blob-num js-line-number" data-line-number="213"></td>
        <td id="LC213" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L214" class="blob-num js-line-number" data-line-number="214"></td>
        <td id="LC214" class="blob-code blob-code-inner js-file-line"><span class="pl-v">calcICL</span> <span class="pl-k">=</span> <span class="pl-k">function</span>(<span class="pl-smi">res</span>,<span class="pl-v">title</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>untitled_consensus_cluster<span class="pl-pds">&quot;</span></span>,<span class="pl-v">plot</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,<span class="pl-v">writeTable</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span>){</td>
      </tr>
      <tr>
        <td id="L215" class="blob-num js-line-number" data-line-number="215"></td>
        <td id="LC215" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>calculates and plots cluster consensus and item consensus</span></td>
      </tr>
      <tr>
        <td id="L216" class="blob-num js-line-number" data-line-number="216"></td>
        <td id="LC216" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">cc</span><span class="pl-k">=</span>rbind()</td>
      </tr>
      <tr>
        <td id="L217" class="blob-num js-line-number" data-line-number="217"></td>
        <td id="LC217" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">cci</span> <span class="pl-k">=</span> rbind()</td>
      </tr>
      <tr>
        <td id="L218" class="blob-num js-line-number" data-line-number="218"></td>
        <td id="LC218" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">sumRes</span><span class="pl-k">=</span><span class="pl-k">list</span>()</td>
      </tr>
      <tr>
        <td id="L219" class="blob-num js-line-number" data-line-number="219"></td>
        <td id="LC219" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">colorsArr</span><span class="pl-k">=</span>c()</td>
      </tr>
      <tr>
        <td id="L220" class="blob-num js-line-number" data-line-number="220"></td>
        <td id="LC220" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L221" class="blob-num js-line-number" data-line-number="221"></td>
        <td id="LC221" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>make results directory</span></td>
      </tr>
      <tr>
        <td id="L222" class="blob-num js-line-number" data-line-number="222"></td>
        <td id="LC222" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>((is.null(<span class="pl-smi">plot</span>)<span class="pl-k">==</span><span class="pl-c1">FALSE</span> <span class="pl-k">|</span> <span class="pl-smi">writeTable</span>) <span class="pl-k">&amp;</span> <span class="pl-k">!</span>file.exists(paste(<span class="pl-smi">title</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))){</td>
      </tr>
      <tr>
        <td id="L223" class="blob-num js-line-number" data-line-number="223"></td>
        <td id="LC223" class="blob-code blob-code-inner js-file-line">	dir.create(paste(<span class="pl-smi">title</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L224" class="blob-num js-line-number" data-line-number="224"></td>
        <td id="LC224" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L225" class="blob-num js-line-number" data-line-number="225"></td>
        <td id="LC225" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">plot</span>)){</td>
      </tr>
      <tr>
        <td id="L226" class="blob-num js-line-number" data-line-number="226"></td>
        <td id="LC226" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>to screen</span></td>
      </tr>
      <tr>
        <td id="L227" class="blob-num js-line-number" data-line-number="227"></td>
        <td id="LC227" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>pdf<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L228" class="blob-num js-line-number" data-line-number="228"></td>
        <td id="LC228" class="blob-code blob-code-inner js-file-line">    pdf(<span class="pl-v">onefile</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>, paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>icl.pdf<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L229" class="blob-num js-line-number" data-line-number="229"></td>
        <td id="LC229" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>ps<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L230" class="blob-num js-line-number" data-line-number="230"></td>
        <td id="LC230" class="blob-code blob-code-inner js-file-line">    postscript(<span class="pl-v">onefile</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>, paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>icl.ps<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L231" class="blob-num js-line-number" data-line-number="231"></td>
        <td id="LC231" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span> (<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>png<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L232" class="blob-num js-line-number" data-line-number="232"></td>
        <td id="LC232" class="blob-code blob-code-inner js-file-line">    png(paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>icl%03d.png<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L233" class="blob-num js-line-number" data-line-number="233"></td>
        <td id="LC233" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span> (<span class="pl-smi">plot</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>pngBMP<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L234" class="blob-num js-line-number" data-line-number="234"></td>
        <td id="LC234" class="blob-code blob-code-inner js-file-line">    bitmap(paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>icl%03d.png<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>))</td>
      </tr>
      <tr>
        <td id="L235" class="blob-num js-line-number" data-line-number="235"></td>
        <td id="LC235" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L236" class="blob-num js-line-number" data-line-number="236"></td>
        <td id="LC236" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L237" class="blob-num js-line-number" data-line-number="237"></td>
        <td id="LC237" class="blob-code blob-code-inner js-file-line">  par(<span class="pl-v">mfrow</span><span class="pl-k">=</span>c(<span class="pl-c1">3</span>,<span class="pl-c1">1</span>),<span class="pl-v">mar</span><span class="pl-k">=</span>c(<span class="pl-c1">4</span>,<span class="pl-c1">3</span>,<span class="pl-c1">2</span>,<span class="pl-c1">0</span>))</td>
      </tr>
      <tr>
        <td id="L238" class="blob-num js-line-number" data-line-number="238"></td>
        <td id="LC238" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L239" class="blob-num js-line-number" data-line-number="239"></td>
        <td id="LC239" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> (<span class="pl-smi">k</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span>length(<span class="pl-smi">res</span>)){ <span class="pl-c"><span class="pl-c">#</span>each k</span></td>
      </tr>
      <tr>
        <td id="L240" class="blob-num js-line-number" data-line-number="240"></td>
        <td id="LC240" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">eiCols</span> <span class="pl-k">=</span> c();</td>
      </tr>
      <tr>
        <td id="L241" class="blob-num js-line-number" data-line-number="241"></td>
        <td id="LC241" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">o</span> <span class="pl-k">=</span> <span class="pl-smi">res</span>[[<span class="pl-smi">k</span>]]</td>
      </tr>
      <tr>
        <td id="L242" class="blob-num js-line-number" data-line-number="242"></td>
        <td id="LC242" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">m</span> <span class="pl-k">=</span> <span class="pl-smi">o</span><span class="pl-k">$</span><span class="pl-smi">consensusMatrix</span></td>
      </tr>
      <tr>
        <td id="L243" class="blob-num js-line-number" data-line-number="243"></td>
        <td id="LC243" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">m</span> <span class="pl-k">=</span> triangle(<span class="pl-smi">m</span>,<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-c1">2</span>)</td>
      </tr>
      <tr>
        <td id="L244" class="blob-num js-line-number" data-line-number="244"></td>
        <td id="LC244" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span> (<span class="pl-smi">ci</span> <span class="pl-k">in</span> sort(unique(<span class="pl-smi">o</span><span class="pl-k">$</span><span class="pl-smi">consensusClass</span>))){ <span class="pl-c"><span class="pl-c">#</span>each cluster in k</span></td>
      </tr>
      <tr>
        <td id="L245" class="blob-num js-line-number" data-line-number="245"></td>
        <td id="LC245" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">items</span> <span class="pl-k">=</span> which(<span class="pl-smi">o</span><span class="pl-k">$</span><span class="pl-smi">consensusClass</span><span class="pl-k">==</span><span class="pl-smi">ci</span>)</td>
      </tr>
      <tr>
        <td id="L246" class="blob-num js-line-number" data-line-number="246"></td>
        <td id="LC246" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">nk</span> <span class="pl-k">=</span> length(<span class="pl-smi">items</span>)</td>
      </tr>
      <tr>
        <td id="L247" class="blob-num js-line-number" data-line-number="247"></td>
        <td id="LC247" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">mk</span> <span class="pl-k">=</span> sum( <span class="pl-smi">m</span>[<span class="pl-smi">items</span>,<span class="pl-smi">items</span>], <span class="pl-v">na.rm</span><span class="pl-k">=</span><span class="pl-c1">T</span>)<span class="pl-k">/</span>((<span class="pl-smi">nk</span><span class="pl-k">*</span>(<span class="pl-smi">nk</span><span class="pl-k">-</span><span class="pl-c1">1</span>))<span class="pl-k">/</span><span class="pl-c1">2</span>)</td>
      </tr>
      <tr>
        <td id="L248" class="blob-num js-line-number" data-line-number="248"></td>
        <td id="LC248" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">cc</span><span class="pl-k">=</span>rbind(<span class="pl-smi">cc</span>,c(<span class="pl-smi">k</span>,<span class="pl-smi">ci</span>,<span class="pl-smi">mk</span>)) <span class="pl-c"><span class="pl-c">#</span>cluster-consensus</span></td>
      </tr>
      <tr>
        <td id="L249" class="blob-num js-line-number" data-line-number="249"></td>
        <td id="LC249" class="blob-code blob-code-inner js-file-line">	</td>
      </tr>
      <tr>
        <td id="L250" class="blob-num js-line-number" data-line-number="250"></td>
        <td id="LC250" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">for</span> (<span class="pl-smi">ei</span> <span class="pl-k">in</span> rev(<span class="pl-smi">res</span>[[<span class="pl-c1">2</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusTree</span><span class="pl-k">$</span><span class="pl-smi">order</span>) ){</td>
      </tr>
      <tr>
        <td id="L251" class="blob-num js-line-number" data-line-number="251"></td>
        <td id="LC251" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">denom</span> <span class="pl-k">=</span> <span class="pl-k">if</span> (<span class="pl-smi">ei</span> <span class="pl-k">%in%</span> <span class="pl-smi">items</span>) { <span class="pl-smi">nk</span> <span class="pl-k">-</span> <span class="pl-c1">1</span>} <span class="pl-k">else</span> { <span class="pl-smi">nk</span> }</td>
      </tr>
      <tr>
        <td id="L252" class="blob-num js-line-number" data-line-number="252"></td>
        <td id="LC252" class="blob-code blob-code-inner js-file-line">        	<span class="pl-v">mei</span> <span class="pl-k">=</span> sum( c(<span class="pl-smi">m</span>[<span class="pl-smi">ei</span>,<span class="pl-smi">items</span>],<span class="pl-smi">m</span>[<span class="pl-smi">items</span>,<span class="pl-smi">ei</span>]), <span class="pl-v">na.rm</span><span class="pl-k">=</span><span class="pl-c1">T</span>)<span class="pl-k">/</span><span class="pl-smi">denom</span>  <span class="pl-c"><span class="pl-c">#</span> mean item consensus to a cluster.</span></td>
      </tr>
      <tr>
        <td id="L253" class="blob-num js-line-number" data-line-number="253"></td>
        <td id="LC253" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">cci</span> <span class="pl-k">=</span> rbind(<span class="pl-smi">cci</span>,c(<span class="pl-smi">k</span>,<span class="pl-smi">ci</span>,<span class="pl-smi">ei</span>,<span class="pl-smi">mei</span>)) <span class="pl-c"><span class="pl-c">#</span>cluster, cluster index, item index, item-consensus</span></td>
      </tr>
      <tr>
        <td id="L254" class="blob-num js-line-number" data-line-number="254"></td>
        <td id="LC254" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L255" class="blob-num js-line-number" data-line-number="255"></td>
        <td id="LC255" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">eiCols</span> <span class="pl-k">=</span> c(<span class="pl-smi">eiCols</span>, rep(<span class="pl-smi">ci</span>,length(<span class="pl-smi">o</span><span class="pl-k">$</span><span class="pl-smi">consensusClass</span>)) )</td>
      </tr>
      <tr>
        <td id="L256" class="blob-num js-line-number" data-line-number="256"></td>
        <td id="LC256" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L257" class="blob-num js-line-number" data-line-number="257"></td>
        <td id="LC257" class="blob-code blob-code-inner js-file-line">	  </td>
      </tr>
      <tr>
        <td id="L258" class="blob-num js-line-number" data-line-number="258"></td>
        <td id="LC258" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">cck</span> <span class="pl-k">=</span> <span class="pl-smi">cci</span>[which(<span class="pl-smi">cci</span>[,<span class="pl-c1">1</span>]<span class="pl-k">==</span><span class="pl-smi">k</span>),] <span class="pl-c"><span class="pl-c">#</span>only plot the new k data.</span></td>
      </tr>
      <tr>
        <td id="L259" class="blob-num js-line-number" data-line-number="259"></td>
        <td id="LC259" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L260" class="blob-num js-line-number" data-line-number="260"></td>
        <td id="LC260" class="blob-code blob-code-inner js-file-line">	  <span class="pl-c"><span class="pl-c">#</span>group by item, order by cluster i</span></td>
      </tr>
      <tr>
        <td id="L261" class="blob-num js-line-number" data-line-number="261"></td>
        <td id="LC261" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">w</span><span class="pl-k">=</span>lapply(split(<span class="pl-smi">cck</span>,<span class="pl-smi">cck</span>[,<span class="pl-c1">3</span>]), <span class="pl-k">function</span>(<span class="pl-smi">x</span>) { <span class="pl-v">y</span><span class="pl-k">=</span><span class="pl-k">matrix</span>(unlist(<span class="pl-smi">x</span>),<span class="pl-v">ncol</span><span class="pl-k">=</span><span class="pl-c1">4</span>); <span class="pl-smi">y</span>[order(<span class="pl-smi">y</span>[,<span class="pl-c1">2</span>]),<span class="pl-c1">4</span>] }) </td>
      </tr>
      <tr>
        <td id="L262" class="blob-num js-line-number" data-line-number="262"></td>
        <td id="LC262" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">q</span> <span class="pl-k">=</span> <span class="pl-k">matrix</span>(as.numeric(unlist(<span class="pl-smi">w</span>)),<span class="pl-v">ncol</span><span class="pl-k">=</span>length(<span class="pl-smi">w</span>),<span class="pl-v">byrow</span><span class="pl-k">=</span><span class="pl-c1">F</span>)</td>
      </tr>
      <tr>
        <td id="L263" class="blob-num js-line-number" data-line-number="263"></td>
        <td id="LC263" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">q</span> <span class="pl-k">=</span> <span class="pl-smi">q</span>[,<span class="pl-smi">res</span>[[<span class="pl-c1">2</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusTree</span><span class="pl-k">$</span><span class="pl-smi">order</span>] <span class="pl-c"><span class="pl-c">#</span>order by leave order of k=2</span></td>
      </tr>
      <tr>
        <td id="L264" class="blob-num js-line-number" data-line-number="264"></td>
        <td id="LC264" class="blob-code blob-code-inner js-file-line"> 	  <span class="pl-c"><span class="pl-c">#</span>q is a matrix of k rows and sample columns, values are item consensus of sample to the cluster.</span></td>
      </tr>
      <tr>
        <td id="L265" class="blob-num js-line-number" data-line-number="265"></td>
        <td id="LC265" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L266" class="blob-num js-line-number" data-line-number="266"></td>
        <td id="LC266" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">thisColors</span> <span class="pl-k">=</span> unique(cbind(<span class="pl-smi">res</span>[[<span class="pl-smi">k</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusClass</span>,<span class="pl-smi">res</span>[[<span class="pl-smi">k</span>]]<span class="pl-k">$</span><span class="pl-smi">clrs</span>[[<span class="pl-c1">1</span>]]))</td>
      </tr>
      <tr>
        <td id="L267" class="blob-num js-line-number" data-line-number="267"></td>
        <td id="LC267" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">thisColors</span><span class="pl-k">=</span><span class="pl-smi">thisColors</span>[order(as.numeric(<span class="pl-smi">thisColors</span>[,<span class="pl-c1">1</span>])),<span class="pl-c1">2</span>]</td>
      </tr>
      <tr>
        <td id="L268" class="blob-num js-line-number" data-line-number="268"></td>
        <td id="LC268" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">colorsArr</span><span class="pl-k">=</span>c(<span class="pl-smi">colorsArr</span>,<span class="pl-smi">thisColors</span>)</td>
      </tr>
      <tr>
        <td id="L269" class="blob-num js-line-number" data-line-number="269"></td>
        <td id="LC269" class="blob-code blob-code-inner js-file-line">	  <span class="pl-smi">sumRes</span>[[<span class="pl-smi">k</span>]] <span class="pl-k">=</span> rankedBarPlot(<span class="pl-smi">q</span>,<span class="pl-smi">thisColors</span>,<span class="pl-v">cc</span><span class="pl-k">=</span><span class="pl-smi">res</span>[[<span class="pl-smi">k</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusClass</span>[<span class="pl-smi">res</span>[[<span class="pl-c1">2</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusTree</span><span class="pl-k">$</span><span class="pl-smi">order</span>],paste(<span class="pl-s"><span class="pl-pds">&quot;</span>k=<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">k</span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>) )</td>
      </tr>
      <tr>
        <td id="L270" class="blob-num js-line-number" data-line-number="270"></td>
        <td id="LC270" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L271" class="blob-num js-line-number" data-line-number="271"></td>
        <td id="LC271" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L272" class="blob-num js-line-number" data-line-number="272"></td>
        <td id="LC272" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">ys</span><span class="pl-k">=</span><span class="pl-v">cs</span><span class="pl-k">=</span><span class="pl-v">lab</span><span class="pl-k">=</span>c()</td>
      </tr>
      <tr>
        <td id="L273" class="blob-num js-line-number" data-line-number="273"></td>
        <td id="LC273" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">lastk</span><span class="pl-k">=</span><span class="pl-smi">cc</span>[<span class="pl-c1">1</span>,<span class="pl-c1">1</span>]</td>
      </tr>
      <tr>
        <td id="L274" class="blob-num js-line-number" data-line-number="274"></td>
        <td id="LC274" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span>(<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span>length(<span class="pl-smi">colorsArr</span>)){</td>
      </tr>
      <tr>
        <td id="L275" class="blob-num js-line-number" data-line-number="275"></td>
        <td id="LC275" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(<span class="pl-smi">lastk</span> <span class="pl-k">!=</span> <span class="pl-smi">cc</span>[<span class="pl-smi">i</span>,<span class="pl-c1">1</span>]){</td>
      </tr>
      <tr>
        <td id="L276" class="blob-num js-line-number" data-line-number="276"></td>
        <td id="LC276" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">ys</span><span class="pl-k">=</span>c(<span class="pl-smi">ys</span>,<span class="pl-c1">0</span>,<span class="pl-c1">0</span>)</td>
      </tr>
      <tr>
        <td id="L277" class="blob-num js-line-number" data-line-number="277"></td>
        <td id="LC277" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">cs</span><span class="pl-k">=</span>c(<span class="pl-smi">cs</span>,<span class="pl-c1">NA</span>,<span class="pl-c1">NA</span>)</td>
      </tr>
      <tr>
        <td id="L278" class="blob-num js-line-number" data-line-number="278"></td>
        <td id="LC278" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">lastk</span><span class="pl-k">=</span><span class="pl-smi">cc</span>[<span class="pl-smi">i</span>,<span class="pl-c1">1</span>]</td>
      </tr>
      <tr>
        <td id="L279" class="blob-num js-line-number" data-line-number="279"></td>
        <td id="LC279" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">lab</span><span class="pl-k">=</span>c(<span class="pl-smi">lab</span>,<span class="pl-c1">NA</span>,<span class="pl-c1">NA</span>)</td>
      </tr>
      <tr>
        <td id="L280" class="blob-num js-line-number" data-line-number="280"></td>
        <td id="LC280" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L281" class="blob-num js-line-number" data-line-number="281"></td>
        <td id="LC281" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">ys</span><span class="pl-k">=</span>c(<span class="pl-smi">ys</span>,<span class="pl-smi">cc</span>[<span class="pl-smi">i</span>,<span class="pl-c1">3</span>])</td>
      </tr>
      <tr>
        <td id="L282" class="blob-num js-line-number" data-line-number="282"></td>
        <td id="LC282" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">cs</span><span class="pl-k">=</span>c(<span class="pl-smi">cs</span>,<span class="pl-smi">colorsArr</span>[<span class="pl-smi">i</span>])</td>
      </tr>
      <tr>
        <td id="L283" class="blob-num js-line-number" data-line-number="283"></td>
        <td id="LC283" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">lab</span><span class="pl-k">=</span>c(<span class="pl-smi">lab</span>,<span class="pl-smi">cc</span>[<span class="pl-smi">i</span>,<span class="pl-c1">1</span>])</td>
      </tr>
      <tr>
        <td id="L284" class="blob-num js-line-number" data-line-number="284"></td>
        <td id="LC284" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L285" class="blob-num js-line-number" data-line-number="285"></td>
        <td id="LC285" class="blob-code blob-code-inner js-file-line">  names(<span class="pl-smi">ys</span>) <span class="pl-k">=</span> <span class="pl-smi">lab</span></td>
      </tr>
      <tr>
        <td id="L286" class="blob-num js-line-number" data-line-number="286"></td>
        <td id="LC286" class="blob-code blob-code-inner js-file-line">  par(<span class="pl-v">mfrow</span><span class="pl-k">=</span>c(<span class="pl-c1">3</span>,<span class="pl-c1">1</span>),<span class="pl-v">mar</span><span class="pl-k">=</span>c(<span class="pl-c1">4</span>,<span class="pl-c1">3</span>,<span class="pl-c1">2</span>,<span class="pl-c1">0</span>))</td>
      </tr>
      <tr>
        <td id="L287" class="blob-num js-line-number" data-line-number="287"></td>
        <td id="LC287" class="blob-code blob-code-inner js-file-line">  barplot(<span class="pl-smi">ys</span>,<span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">cs</span>,<span class="pl-v">border</span><span class="pl-k">=</span><span class="pl-smi">cs</span>,<span class="pl-v">main</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>cluster-consensus<span class="pl-pds">&quot;</span></span>,<span class="pl-v">ylim</span><span class="pl-k">=</span>c(<span class="pl-c1">0</span>,<span class="pl-c1">1</span>),<span class="pl-v">las</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L288" class="blob-num js-line-number" data-line-number="288"></td>
        <td id="LC288" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(is.null(<span class="pl-smi">plot</span>)<span class="pl-k">==</span><span class="pl-c1">FALSE</span>){</td>
      </tr>
      <tr>
        <td id="L289" class="blob-num js-line-number" data-line-number="289"></td>
        <td id="LC289" class="blob-code blob-code-inner js-file-line">	  dev.off()</td>
      </tr>
      <tr>
        <td id="L290" class="blob-num js-line-number" data-line-number="290"></td>
        <td id="LC290" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L291" class="blob-num js-line-number" data-line-number="291"></td>
        <td id="LC291" class="blob-code blob-code-inner js-file-line">  colnames(<span class="pl-smi">cc</span>) <span class="pl-k">=</span> c(<span class="pl-s"><span class="pl-pds">&quot;</span>k<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>cluster<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>clusterConsensus<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L292" class="blob-num js-line-number" data-line-number="292"></td>
        <td id="LC292" class="blob-code blob-code-inner js-file-line">  colnames(<span class="pl-smi">cci</span>) <span class="pl-k">=</span> c(<span class="pl-s"><span class="pl-pds">&quot;</span>k<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>cluster<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>item<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>itemConsensus<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L293" class="blob-num js-line-number" data-line-number="293"></td>
        <td id="LC293" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">cci</span>[,<span class="pl-s"><span class="pl-pds">&quot;</span>item<span class="pl-pds">&quot;</span></span>] <span class="pl-k">=</span> names(<span class="pl-smi">res</span>[[<span class="pl-c1">2</span>]]<span class="pl-k">$</span><span class="pl-smi">consensusClass</span>)[ <span class="pl-smi">cci</span>[,<span class="pl-s"><span class="pl-pds">&quot;</span>item<span class="pl-pds">&quot;</span></span>] ]</td>
      </tr>
      <tr>
        <td id="L294" class="blob-num js-line-number" data-line-number="294"></td>
        <td id="LC294" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>type cci</span></td>
      </tr>
      <tr>
        <td id="L295" class="blob-num js-line-number" data-line-number="295"></td>
        <td id="LC295" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">cci</span> <span class="pl-k">=</span> <span class="pl-k">data.frame</span>( <span class="pl-v">k</span><span class="pl-k">=</span>as.numeric(<span class="pl-smi">cci</span>[,<span class="pl-s"><span class="pl-pds">&quot;</span>k<span class="pl-pds">&quot;</span></span>]), <span class="pl-v">cluster</span><span class="pl-k">=</span>as.numeric(<span class="pl-smi">cci</span>[,<span class="pl-s"><span class="pl-pds">&quot;</span>cluster<span class="pl-pds">&quot;</span></span>]), <span class="pl-v">item</span><span class="pl-k">=</span><span class="pl-smi">cci</span>[,<span class="pl-s"><span class="pl-pds">&quot;</span>item<span class="pl-pds">&quot;</span></span>], <span class="pl-v">itemConsensus</span><span class="pl-k">=</span>as.numeric(<span class="pl-smi">cci</span>[,<span class="pl-s"><span class="pl-pds">&quot;</span>itemConsensus<span class="pl-pds">&quot;</span></span>])) </td>
      </tr>
      <tr>
        <td id="L296" class="blob-num js-line-number" data-line-number="296"></td>
        <td id="LC296" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L297" class="blob-num js-line-number" data-line-number="297"></td>
        <td id="LC297" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>write to file.</span></td>
      </tr>
      <tr>
        <td id="L298" class="blob-num js-line-number" data-line-number="298"></td>
        <td id="LC298" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(<span class="pl-smi">writeTable</span>){</td>
      </tr>
      <tr>
        <td id="L299" class="blob-num js-line-number" data-line-number="299"></td>
        <td id="LC299" class="blob-code blob-code-inner js-file-line">	write.csv(<span class="pl-v">file</span><span class="pl-k">=</span>paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.summary.cluster.consensus.csv<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>),<span class="pl-v">row.names</span><span class="pl-k">=</span><span class="pl-c1">F</span>, <span class="pl-smi">cc</span>)</td>
      </tr>
      <tr>
        <td id="L300" class="blob-num js-line-number" data-line-number="300"></td>
        <td id="LC300" class="blob-code blob-code-inner js-file-line">	write.csv(<span class="pl-v">file</span><span class="pl-k">=</span>paste(<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>/<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">title</span>,<span class="pl-s"><span class="pl-pds">&quot;</span>.summary.item.consensus.csv<span class="pl-pds">&quot;</span></span>,<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>), <span class="pl-v">row.names</span><span class="pl-k">=</span><span class="pl-c1">F</span>, <span class="pl-smi">cc</span>)</td>
      </tr>
      <tr>
        <td id="L301" class="blob-num js-line-number" data-line-number="301"></td>
        <td id="LC301" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L302" class="blob-num js-line-number" data-line-number="302"></td>
        <td id="LC302" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">return</span>(<span class="pl-k">list</span>(<span class="pl-v">clusterConsensus</span><span class="pl-k">=</span><span class="pl-smi">cc</span>,<span class="pl-v">itemConsensus</span><span class="pl-k">=</span><span class="pl-smi">cci</span>))</td>
      </tr>
      <tr>
        <td id="L303" class="blob-num js-line-number" data-line-number="303"></td>
        <td id="LC303" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L304" class="blob-num js-line-number" data-line-number="304"></td>
        <td id="LC304" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L305" class="blob-num js-line-number" data-line-number="305"></td>
        <td id="LC305" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L306" class="blob-num js-line-number" data-line-number="306"></td>
        <td id="LC306" class="blob-code blob-code-inner js-file-line"><span class="pl-en">ccRun</span> <span class="pl-k">&lt;-</span> <span class="pl-k">function</span>( <span class="pl-v">d</span><span class="pl-k">=</span><span class="pl-smi">d</span>,</td>
      </tr>
      <tr>
        <td id="L307" class="blob-num js-line-number" data-line-number="307"></td>
        <td id="LC307" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">maxK</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L308" class="blob-num js-line-number" data-line-number="308"></td>
        <td id="LC308" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">repCount</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L309" class="blob-num js-line-number" data-line-number="309"></td>
        <td id="LC309" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">diss</span><span class="pl-k">=</span>inherits( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>dist<span class="pl-pds">&quot;</span></span> ),</td>
      </tr>
      <tr>
        <td id="L310" class="blob-num js-line-number" data-line-number="310"></td>
        <td id="LC310" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">pItem</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L311" class="blob-num js-line-number" data-line-number="311"></td>
        <td id="LC311" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">pFeature</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L312" class="blob-num js-line-number" data-line-number="312"></td>
        <td id="LC312" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">innerLinkage</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L313" class="blob-num js-line-number" data-line-number="313"></td>
        <td id="LC313" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">distance</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>, <span class="pl-c"><span class="pl-c">#</span>ifelse( inherits(d,&quot;dist&quot;), attr( d, &quot;method&quot; ), &quot;euclidean&quot; ),@@@@@</span></td>
      </tr>
      <tr>
        <td id="L314" class="blob-num js-line-number" data-line-number="314"></td>
        <td id="LC314" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">clusterAlg</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L315" class="blob-num js-line-number" data-line-number="315"></td>
        <td id="LC315" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">weightsItem</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L316" class="blob-num js-line-number" data-line-number="316"></td>
        <td id="LC316" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">weightsFeature</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L317" class="blob-num js-line-number" data-line-number="317"></td>
        <td id="LC317" class="blob-code blob-code-inner js-file-line">                   <span class="pl-v">verbose</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L318" class="blob-num js-line-number" data-line-number="318"></td>
        <td id="LC318" class="blob-code blob-code-inner js-file-line">		   <span class="pl-v">corUse</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>) {</td>
      </tr>
      <tr>
        <td id="L319" class="blob-num js-line-number" data-line-number="319"></td>
        <td id="LC319" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">m</span> <span class="pl-k">=</span> vector(<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>list<span class="pl-pds">&#39;</span></span>, <span class="pl-smi">repCount</span>)</td>
      </tr>
      <tr>
        <td id="L320" class="blob-num js-line-number" data-line-number="320"></td>
        <td id="LC320" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">ml</span> <span class="pl-k">=</span> vector(<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>list<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">maxK</span>)</td>
      </tr>
      <tr>
        <td id="L321" class="blob-num js-line-number" data-line-number="321"></td>
        <td id="LC321" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">n</span> <span class="pl-k">&lt;-</span> ifelse( <span class="pl-smi">diss</span>, ncol( as.matrix(<span class="pl-smi">d</span>) ), ncol(<span class="pl-smi">d</span>) )</td>
      </tr>
      <tr>
        <td id="L322" class="blob-num js-line-number" data-line-number="322"></td>
        <td id="LC322" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">mCount</span> <span class="pl-k">=</span> <span class="pl-v">mConsist</span> <span class="pl-k">=</span> <span class="pl-k">matrix</span>(c(<span class="pl-c1">0</span>),<span class="pl-v">ncol</span><span class="pl-k">=</span><span class="pl-smi">n</span>,<span class="pl-v">nrow</span><span class="pl-k">=</span><span class="pl-smi">n</span>)</td>
      </tr>
      <tr>
        <td id="L323" class="blob-num js-line-number" data-line-number="323"></td>
        <td id="LC323" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">ml</span>[[<span class="pl-c1">1</span>]] <span class="pl-k">=</span> c(<span class="pl-c1">0</span>);</td>
      </tr>
      <tr>
        <td id="L324" class="blob-num js-line-number" data-line-number="324"></td>
        <td id="LC324" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L325" class="blob-num js-line-number" data-line-number="325"></td>
        <td id="LC325" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span> (is.null( <span class="pl-smi">distance</span> ) ) <span class="pl-smi">distance</span> <span class="pl-k">&lt;-</span> <span class="pl-s"><span class="pl-pds">&#39;</span>euclidean<span class="pl-pds">&#39;</span></span>  <span class="pl-c"><span class="pl-c">#</span># necessary if d is a dist object and attr( d, &quot;method&quot; ) == NULL</span></td>
      </tr>
      <tr>
        <td id="L326" class="blob-num js-line-number" data-line-number="326"></td>
        <td id="LC326" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L327" class="blob-num js-line-number" data-line-number="327"></td>
        <td id="LC327" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">acceptable.distance</span> <span class="pl-k">&lt;-</span> c( <span class="pl-s"><span class="pl-pds">&quot;</span>euclidean<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>maximum<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>manhattan<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>canberra<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>binary<span class="pl-pds">&quot;</span></span>,<span class="pl-s"><span class="pl-pds">&quot;</span>minkowski<span class="pl-pds">&quot;</span></span>,</td>
      </tr>
      <tr>
        <td id="L328" class="blob-num js-line-number" data-line-number="328"></td>
        <td id="LC328" class="blob-code blob-code-inner js-file-line">                      <span class="pl-s"><span class="pl-pds">&quot;</span>pearson<span class="pl-pds">&quot;</span></span>, <span class="pl-s"><span class="pl-pds">&quot;</span>spearman<span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L329" class="blob-num js-line-number" data-line-number="329"></td>
        <td id="LC329" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L330" class="blob-num js-line-number" data-line-number="330"></td>
        <td id="LC330" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">main.dist.obj</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">NULL</span></td>
      </tr>
      <tr>
        <td id="L331" class="blob-num js-line-number" data-line-number="331"></td>
        <td id="LC331" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span> ( <span class="pl-smi">diss</span> ){</td>
      </tr>
      <tr>
        <td id="L332" class="blob-num js-line-number" data-line-number="332"></td>
        <td id="LC332" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">main.dist.obj</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">d</span></td>
      </tr>
      <tr>
        <td id="L333" class="blob-num js-line-number" data-line-number="333"></td>
        <td id="LC333" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span># reset the pFeature &amp; weightsFeature params if they&#39;ve been set (irrelevant if d is a dist matrix)</span></td>
      </tr>
      <tr>
        <td id="L334" class="blob-num js-line-number" data-line-number="334"></td>
        <td id="LC334" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( ( <span class="pl-k">!</span>is.null(<span class="pl-smi">pFeature</span>) ) <span class="pl-k">&amp;&amp;</span></td>
      </tr>
      <tr>
        <td id="L335" class="blob-num js-line-number" data-line-number="335"></td>
        <td id="LC335" class="blob-code blob-code-inner js-file-line">         ( <span class="pl-smi">pFeature</span> <span class="pl-k">&lt;</span> <span class="pl-c1">1</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L336" class="blob-num js-line-number" data-line-number="336"></td>
        <td id="LC336" class="blob-code blob-code-inner js-file-line">      message( <span class="pl-s"><span class="pl-pds">&quot;</span>user-supplied data is a distance matrix; ignoring user-specified pFeature parameter<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L337" class="blob-num js-line-number" data-line-number="337"></td>
        <td id="LC337" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">pFeature</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">1</span> <span class="pl-c"><span class="pl-c">#</span> set it to 1 to avoid problems with sampleCols</span></td>
      </tr>
      <tr>
        <td id="L338" class="blob-num js-line-number" data-line-number="338"></td>
        <td id="LC338" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L339" class="blob-num js-line-number" data-line-number="339"></td>
        <td id="LC339" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( <span class="pl-k">!</span> is.null( <span class="pl-smi">weightsFeature</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L340" class="blob-num js-line-number" data-line-number="340"></td>
        <td id="LC340" class="blob-code blob-code-inner js-file-line">      message( <span class="pl-s"><span class="pl-pds">&quot;</span>user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter<span class="pl-cce">\n</span><span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L341" class="blob-num js-line-number" data-line-number="341"></td>
        <td id="LC341" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">weightsFeature</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">NULL</span>  <span class="pl-c"><span class="pl-c">#</span> set it to NULL to avoid problems with sampleCols</span></td>
      </tr>
      <tr>
        <td id="L342" class="blob-num js-line-number" data-line-number="342"></td>
        <td id="LC342" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L343" class="blob-num js-line-number" data-line-number="343"></td>
        <td id="LC343" class="blob-code blob-code-inner js-file-line">  } <span class="pl-k">else</span> { <span class="pl-c"><span class="pl-c">#</span># d is a data matrix</span></td>
      </tr>
      <tr>
        <td id="L344" class="blob-num js-line-number" data-line-number="344"></td>
        <td id="LC344" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span># we&#39;re not sampling over the features</span></td>
      </tr>
      <tr>
        <td id="L345" class="blob-num js-line-number" data-line-number="345"></td>
        <td id="LC345" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( ( <span class="pl-smi">clusterAlg</span> <span class="pl-k">!=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>km<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&amp;&amp;</span></td>
      </tr>
      <tr>
        <td id="L346" class="blob-num js-line-number" data-line-number="346"></td>
        <td id="LC346" class="blob-code blob-code-inner js-file-line">         ( is.null( <span class="pl-smi">pFeature</span> ) <span class="pl-k">||</span></td>
      </tr>
      <tr>
        <td id="L347" class="blob-num js-line-number" data-line-number="347"></td>
        <td id="LC347" class="blob-code blob-code-inner js-file-line">           ( ( <span class="pl-smi">pFeature</span> <span class="pl-k">==</span> <span class="pl-c1">1</span> ) <span class="pl-k">&amp;&amp;</span> is.null( <span class="pl-smi">weightsFeature</span> ) ) ) ) {</td>
      </tr>
      <tr>
        <td id="L348" class="blob-num js-line-number" data-line-number="348"></td>
        <td id="LC348" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># only generate a main.dist.object IFF 1) d is a matrix, 2) we&#39;re not sampling the features, and 3) the algorithm isn&#39;t &#39;km&#39;</span></td>
      </tr>
      <tr>
        <td id="L349" class="blob-num js-line-number" data-line-number="349"></td>
        <td id="LC349" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( inherits( <span class="pl-smi">distance</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>character<span class="pl-pds">&quot;</span></span> ) ) {</td>
      </tr>
      <tr>
        <td id="L350" class="blob-num js-line-number" data-line-number="350"></td>
        <td id="LC350" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">if</span> ( <span class="pl-k">!</span> <span class="pl-smi">distance</span> <span class="pl-k">%in%</span>  <span class="pl-smi">acceptable.distance</span>  <span class="pl-k">&amp;</span>  ( class(try(get(<span class="pl-smi">distance</span>),<span class="pl-v">silent</span><span class="pl-k">=</span><span class="pl-c1">T</span>))<span class="pl-k">!=</span><span class="pl-s"><span class="pl-pds">&quot;</span>function<span class="pl-pds">&quot;</span></span>) ) stop(<span class="pl-s"><span class="pl-pds">&quot;</span>unsupported distance.<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L351" class="blob-num js-line-number" data-line-number="351"></td>
        <td id="LC351" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L352" class="blob-num js-line-number" data-line-number="352"></td>
        <td id="LC352" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">if</span>(<span class="pl-smi">distance</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>pearson<span class="pl-pds">&quot;</span></span> <span class="pl-k">|</span> <span class="pl-smi">distance</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>spearman<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L353" class="blob-num js-line-number" data-line-number="353"></td>
        <td id="LC353" class="blob-code blob-code-inner js-file-line">          <span class="pl-smi">main.dist.obj</span> <span class="pl-k">&lt;-</span> as.dist( <span class="pl-c1">1</span><span class="pl-k">-</span>cor(<span class="pl-smi">d</span>,<span class="pl-v">method</span><span class="pl-k">=</span><span class="pl-smi">distance</span>,<span class="pl-v">use</span><span class="pl-k">=</span><span class="pl-smi">corUse</span> ))</td>
      </tr>
      <tr>
        <td id="L354" class="blob-num js-line-number" data-line-number="354"></td>
        <td id="LC354" class="blob-code blob-code-inner js-file-line">        }<span class="pl-k">else</span> <span class="pl-k">if</span>( class(try(get(<span class="pl-smi">distance</span>),<span class="pl-v">silent</span><span class="pl-k">=</span><span class="pl-c1">T</span>))<span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>function<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L355" class="blob-num js-line-number" data-line-number="355"></td>
        <td id="LC355" class="blob-code blob-code-inner js-file-line">          <span class="pl-smi">main.dist.obj</span> <span class="pl-k">&lt;-</span> get(<span class="pl-smi">distance</span>)( t( <span class="pl-smi">d</span> )   )</td>
      </tr>
      <tr>
        <td id="L356" class="blob-num js-line-number" data-line-number="356"></td>
        <td id="LC356" class="blob-code blob-code-inner js-file-line">	}<span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L357" class="blob-num js-line-number" data-line-number="357"></td>
        <td id="LC357" class="blob-code blob-code-inner js-file-line">          <span class="pl-smi">main.dist.obj</span> <span class="pl-k">&lt;-</span> dist( t(<span class="pl-smi">d</span>), <span class="pl-v">method</span><span class="pl-k">=</span><span class="pl-smi">distance</span> )</td>
      </tr>
      <tr>
        <td id="L358" class="blob-num js-line-number" data-line-number="358"></td>
        <td id="LC358" class="blob-code blob-code-inner js-file-line">        }</td>
      </tr>
      <tr>
        <td id="L359" class="blob-num js-line-number" data-line-number="359"></td>
        <td id="LC359" class="blob-code blob-code-inner js-file-line">        attr( <span class="pl-smi">main.dist.obj</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&lt;-</span> <span class="pl-smi">distance</span>  </td>
      </tr>
      <tr>
        <td id="L360" class="blob-num js-line-number" data-line-number="360"></td>
        <td id="LC360" class="blob-code blob-code-inner js-file-line">      } <span class="pl-k">else</span> stop(<span class="pl-s"><span class="pl-pds">&quot;</span>unsupported distance specified.<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L361" class="blob-num js-line-number" data-line-number="361"></td>
        <td id="LC361" class="blob-code blob-code-inner js-file-line">    } <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L362" class="blob-num js-line-number" data-line-number="362"></td>
        <td id="LC362" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># pFeature &lt; 1 or a weightsFeature != NULL</span></td>
      </tr>
      <tr>
        <td id="L363" class="blob-num js-line-number" data-line-number="363"></td>
        <td id="LC363" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL</span></td>
      </tr>
      <tr>
        <td id="L364" class="blob-num js-line-number" data-line-number="364"></td>
        <td id="LC364" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L365" class="blob-num js-line-number" data-line-number="365"></td>
        <td id="LC365" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L366" class="blob-num js-line-number" data-line-number="366"></td>
        <td id="LC366" class="blob-code blob-code-inner js-file-line"> </td>
      </tr>
      <tr>
        <td id="L367" class="blob-num js-line-number" data-line-number="367"></td>
        <td id="LC367" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L368" class="blob-num js-line-number" data-line-number="368"></td>
        <td id="LC368" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> (<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span><span class="pl-smi">repCount</span>){</td>
      </tr>
      <tr>
        <td id="L369" class="blob-num js-line-number" data-line-number="369"></td>
        <td id="LC369" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(<span class="pl-smi">verbose</span>){</td>
      </tr>
      <tr>
        <td id="L370" class="blob-num js-line-number" data-line-number="370"></td>
        <td id="LC370" class="blob-code blob-code-inner js-file-line">      message(paste(<span class="pl-s"><span class="pl-pds">&quot;</span>random subsample<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">i</span>));</td>
      </tr>
      <tr>
        <td id="L371" class="blob-num js-line-number" data-line-number="371"></td>
        <td id="LC371" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L372" class="blob-num js-line-number" data-line-number="372"></td>
        <td id="LC372" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span># take expression matrix sample, samples and genes</span></td>
      </tr>
      <tr>
        <td id="L373" class="blob-num js-line-number" data-line-number="373"></td>
        <td id="LC373" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">sample_x</span> <span class="pl-k">=</span> sampleCols( <span class="pl-smi">d</span>, <span class="pl-smi">pItem</span>, <span class="pl-smi">pFeature</span>, <span class="pl-smi">weightsItem</span>, <span class="pl-smi">weightsFeature</span> )</td>
      </tr>
      <tr>
        <td id="L374" class="blob-num js-line-number" data-line-number="374"></td>
        <td id="LC374" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L375" class="blob-num js-line-number" data-line-number="375"></td>
        <td id="LC375" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">this_dist</span> <span class="pl-k">=</span> <span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L376" class="blob-num js-line-number" data-line-number="376"></td>
        <td id="LC376" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( <span class="pl-k">!</span> is.null( <span class="pl-smi">main.dist.obj</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L377" class="blob-num js-line-number" data-line-number="377"></td>
        <td id="LC377" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">boot.cols</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">subcols</span></td>
      </tr>
      <tr>
        <td id="L378" class="blob-num js-line-number" data-line-number="378"></td>
        <td id="LC378" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> as.matrix( <span class="pl-smi">main.dist.obj</span> )[ <span class="pl-smi">boot.cols</span>, <span class="pl-smi">boot.cols</span> ]</td>
      </tr>
      <tr>
        <td id="L379" class="blob-num js-line-number" data-line-number="379"></td>
        <td id="LC379" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( <span class="pl-smi">clusterAlg</span> <span class="pl-k">!=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>km<span class="pl-pds">&quot;</span></span> ) {</td>
      </tr>
      <tr>
        <td id="L380" class="blob-num js-line-number" data-line-number="380"></td>
        <td id="LC380" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span># if this isn&#39;t kmeans, then convert to a distance object</span></td>
      </tr>
      <tr>
        <td id="L381" class="blob-num js-line-number" data-line-number="381"></td>
        <td id="LC381" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> as.dist( <span class="pl-smi">this_dist</span> )</td>
      </tr>
      <tr>
        <td id="L382" class="blob-num js-line-number" data-line-number="382"></td>
        <td id="LC382" class="blob-code blob-code-inner js-file-line">        attr( <span class="pl-smi">this_dist</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&lt;-</span> attr( <span class="pl-smi">main.dist.obj</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L383" class="blob-num js-line-number" data-line-number="383"></td>
        <td id="LC383" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L384" class="blob-num js-line-number" data-line-number="384"></td>
        <td id="LC384" class="blob-code blob-code-inner js-file-line">    } <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L385" class="blob-num js-line-number" data-line-number="385"></td>
        <td id="LC385" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># if main.dist.obj is NULL, then d is a data matrix, and either:</span></td>
      </tr>
      <tr>
        <td id="L386" class="blob-num js-line-number" data-line-number="386"></td>
        <td id="LC386" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#   1) clusterAlg is &#39;km&#39;</span></td>
      </tr>
      <tr>
        <td id="L387" class="blob-num js-line-number" data-line-number="387"></td>
        <td id="LC387" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#   2) pFeatures &lt; 1 or weightsFeatures have been specified, or</span></td>
      </tr>
      <tr>
        <td id="L388" class="blob-num js-line-number" data-line-number="388"></td>
        <td id="LC388" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#   3) both</span></td>
      </tr>
      <tr>
        <td id="L389" class="blob-num js-line-number" data-line-number="389"></td>
        <td id="LC389" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># so we can&#39;t use a main distance object and for every iteration, we will have to re-calculate either</span></td>
      </tr>
      <tr>
        <td id="L390" class="blob-num js-line-number" data-line-number="390"></td>
        <td id="LC390" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#   1) the distance matrix (because we&#39;re also sampling the features as well), or</span></td>
      </tr>
      <tr>
        <td id="L391" class="blob-num js-line-number" data-line-number="391"></td>
        <td id="LC391" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#   2) the submat (if using km) </span></td>
      </tr>
      <tr>
        <td id="L392" class="blob-num js-line-number" data-line-number="392"></td>
        <td id="LC392" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L393" class="blob-num js-line-number" data-line-number="393"></td>
        <td id="LC393" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> ( <span class="pl-smi">clusterAlg</span> <span class="pl-k">!=</span> <span class="pl-s"><span class="pl-pds">&quot;</span>km<span class="pl-pds">&quot;</span></span> )  {</td>
      </tr>
      <tr>
        <td id="L394" class="blob-num js-line-number" data-line-number="394"></td>
        <td id="LC394" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">if</span> ( <span class="pl-k">!</span> <span class="pl-smi">distance</span> <span class="pl-k">%in%</span> <span class="pl-smi">acceptable.distance</span> <span class="pl-k">&amp;</span>  ( class(try(get(<span class="pl-smi">distance</span>),<span class="pl-v">silent</span><span class="pl-k">=</span><span class="pl-c1">T</span>))<span class="pl-k">!=</span><span class="pl-s"><span class="pl-pds">&quot;</span>function<span class="pl-pds">&quot;</span></span>)  ) stop(<span class="pl-s"><span class="pl-pds">&quot;</span>unsupported distance.<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L395" class="blob-num js-line-number" data-line-number="395"></td>
        <td id="LC395" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">if</span>( ( class(try(get(<span class="pl-smi">distance</span>),<span class="pl-v">silent</span><span class="pl-k">=</span><span class="pl-c1">T</span>))<span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>function<span class="pl-pds">&quot;</span></span>) ){</td>
      </tr>
      <tr>
        <td id="L396" class="blob-num js-line-number" data-line-number="396"></td>
        <td id="LC396" class="blob-code blob-code-inner js-file-line">          <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> get(<span class="pl-smi">distance</span>)( t( <span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">submat</span> ) )</td>
      </tr>
      <tr>
        <td id="L397" class="blob-num js-line-number" data-line-number="397"></td>
        <td id="LC397" class="blob-code blob-code-inner js-file-line">	}<span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L398" class="blob-num js-line-number" data-line-number="398"></td>
        <td id="LC398" class="blob-code blob-code-inner js-file-line">	  <span class="pl-k">if</span>( <span class="pl-smi">distance</span> <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>pearson<span class="pl-pds">&quot;</span></span> <span class="pl-k">|</span> <span class="pl-smi">distance</span> <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>spearman<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L399" class="blob-num js-line-number" data-line-number="399"></td>
        <td id="LC399" class="blob-code blob-code-inner js-file-line">            <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> as.dist( <span class="pl-c1">1</span><span class="pl-k">-</span>cor(<span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">submat</span>,<span class="pl-v">use</span><span class="pl-k">=</span><span class="pl-smi">corUse</span>,<span class="pl-v">method</span><span class="pl-k">=</span><span class="pl-smi">distance</span>) )</td>
      </tr>
      <tr>
        <td id="L400" class="blob-num js-line-number" data-line-number="400"></td>
        <td id="LC400" class="blob-code blob-code-inner js-file-line">	  }<span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L401" class="blob-num js-line-number" data-line-number="401"></td>
        <td id="LC401" class="blob-code blob-code-inner js-file-line">            <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> dist( t( <span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">submat</span> ), <span class="pl-v">method</span><span class="pl-k">=</span> <span class="pl-smi">distance</span>  )</td>
      </tr>
      <tr>
        <td id="L402" class="blob-num js-line-number" data-line-number="402"></td>
        <td id="LC402" class="blob-code blob-code-inner js-file-line">	  }</td>
      </tr>
      <tr>
        <td id="L403" class="blob-num js-line-number" data-line-number="403"></td>
        <td id="LC403" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L404" class="blob-num js-line-number" data-line-number="404"></td>
        <td id="LC404" class="blob-code blob-code-inner js-file-line">        attr( <span class="pl-smi">this_dist</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>method<span class="pl-pds">&quot;</span></span> ) <span class="pl-k">&lt;-</span> <span class="pl-smi">distance</span>  </td>
      </tr>
      <tr>
        <td id="L405" class="blob-num js-line-number" data-line-number="405"></td>
        <td id="LC405" class="blob-code blob-code-inner js-file-line">      } <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L406" class="blob-num js-line-number" data-line-number="406"></td>
        <td id="LC406" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span># if we&#39;re not sampling the features, then grab the colslice</span></td>
      </tr>
      <tr>
        <td id="L407" class="blob-num js-line-number" data-line-number="407"></td>
        <td id="LC407" class="blob-code blob-code-inner js-file-line">        <span class="pl-k">if</span> ( is.null( <span class="pl-smi">pFeature</span> ) <span class="pl-k">||</span></td>
      </tr>
      <tr>
        <td id="L408" class="blob-num js-line-number" data-line-number="408"></td>
        <td id="LC408" class="blob-code blob-code-inner js-file-line">            ( ( <span class="pl-smi">pFeature</span> <span class="pl-k">==</span> <span class="pl-c1">1</span> ) <span class="pl-k">&amp;&amp;</span> is.null( <span class="pl-smi">weightsFeature</span> ) ) ) {</td>
      </tr>
      <tr>
        <td id="L409" class="blob-num js-line-number" data-line-number="409"></td>
        <td id="LC409" class="blob-code blob-code-inner js-file-line">          <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">d</span>[, <span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">subcols</span> ]</td>
      </tr>
      <tr>
        <td id="L410" class="blob-num js-line-number" data-line-number="410"></td>
        <td id="LC410" class="blob-code blob-code-inner js-file-line">        } <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L411" class="blob-num js-line-number" data-line-number="411"></td>
        <td id="LC411" class="blob-code blob-code-inner js-file-line">          <span class="pl-k">if</span> ( is.na( <span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">submat</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L412" class="blob-num js-line-number" data-line-number="412"></td>
        <td id="LC412" class="blob-code blob-code-inner js-file-line">            stop( <span class="pl-s"><span class="pl-pds">&quot;</span>error submat is NA<span class="pl-pds">&quot;</span></span> )</td>
      </tr>
      <tr>
        <td id="L413" class="blob-num js-line-number" data-line-number="413"></td>
        <td id="LC413" class="blob-code blob-code-inner js-file-line">          }</td>
      </tr>
      <tr>
        <td id="L414" class="blob-num js-line-number" data-line-number="414"></td>
        <td id="LC414" class="blob-code blob-code-inner js-file-line">          </td>
      </tr>
      <tr>
        <td id="L415" class="blob-num js-line-number" data-line-number="415"></td>
        <td id="LC415" class="blob-code blob-code-inner js-file-line">          <span class="pl-smi">this_dist</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">sample_x</span><span class="pl-k">$</span><span class="pl-smi">submat</span></td>
      </tr>
      <tr>
        <td id="L416" class="blob-num js-line-number" data-line-number="416"></td>
        <td id="LC416" class="blob-code blob-code-inner js-file-line">        } </td>
      </tr>
      <tr>
        <td id="L417" class="blob-num js-line-number" data-line-number="417"></td>
        <td id="LC417" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L418" class="blob-num js-line-number" data-line-number="418"></td>
        <td id="LC418" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L419" class="blob-num js-line-number" data-line-number="419"></td>
        <td id="LC419" class="blob-code blob-code-inner js-file-line">                  </td>
      </tr>
      <tr>
        <td id="L420" class="blob-num js-line-number" data-line-number="420"></td>
        <td id="LC420" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span># cluster samples for HC.</span></td>
      </tr>
      <tr>
        <td id="L421" class="blob-num js-line-number" data-line-number="421"></td>
        <td id="LC421" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">this_cluster</span><span class="pl-k">=</span><span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L422" class="blob-num js-line-number" data-line-number="422"></td>
        <td id="LC422" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span>(<span class="pl-smi">clusterAlg</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>hc<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L423" class="blob-num js-line-number" data-line-number="423"></td>
        <td id="LC423" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">this_cluster</span> <span class="pl-k">=</span> hclust( <span class="pl-smi">this_dist</span>, <span class="pl-v">method</span><span class="pl-k">=</span><span class="pl-smi">innerLinkage</span>)</td>
      </tr>
      <tr>
        <td id="L424" class="blob-num js-line-number" data-line-number="424"></td>
        <td id="LC424" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L425" class="blob-num js-line-number" data-line-number="425"></td>
        <td id="LC425" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>#mCount is possible number of times that two sample occur in same random sample, independent of k</span></td>
      </tr>
      <tr>
        <td id="L426" class="blob-num js-line-number" data-line-number="426"></td>
        <td id="LC426" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>#mCount stores number of times a sample pair was sampled together.</span></td>
      </tr>
      <tr>
        <td id="L427" class="blob-num js-line-number" data-line-number="427"></td>
        <td id="LC427" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">mCount</span> <span class="pl-k">&lt;-</span> connectivityMatrix( rep( <span class="pl-c1">1</span>,length(<span class="pl-smi">sample_x</span>[[<span class="pl-c1">3</span>]])),</td>
      </tr>
      <tr>
        <td id="L428" class="blob-num js-line-number" data-line-number="428"></td>
        <td id="LC428" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-smi">mCount</span>,</td>
      </tr>
      <tr>
        <td id="L429" class="blob-num js-line-number" data-line-number="429"></td>
        <td id="LC429" class="blob-code blob-code-inner js-file-line">                                  <span class="pl-smi">sample_x</span>[[<span class="pl-c1">3</span>]] ) </td>
      </tr>
      <tr>
        <td id="L430" class="blob-num js-line-number" data-line-number="430"></td>
        <td id="LC430" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L431" class="blob-num js-line-number" data-line-number="431"></td>
        <td id="LC431" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>#use samples for each k		</span></td>
      </tr>
      <tr>
        <td id="L432" class="blob-num js-line-number" data-line-number="432"></td>
        <td id="LC432" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span> (<span class="pl-smi">k</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span><span class="pl-smi">maxK</span>){</td>
      </tr>
      <tr>
        <td id="L433" class="blob-num js-line-number" data-line-number="433"></td>
        <td id="LC433" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span>(<span class="pl-smi">verbose</span>){</td>
      </tr>
      <tr>
        <td id="L434" class="blob-num js-line-number" data-line-number="434"></td>
        <td id="LC434" class="blob-code blob-code-inner js-file-line">        message(paste(<span class="pl-s"><span class="pl-pds">&quot;</span>  k =<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">k</span>))</td>
      </tr>
      <tr>
        <td id="L435" class="blob-num js-line-number" data-line-number="435"></td>
        <td id="LC435" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L436" class="blob-num js-line-number" data-line-number="436"></td>
        <td id="LC436" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span> (<span class="pl-smi">i</span><span class="pl-k">==</span><span class="pl-c1">1</span>){</td>
      </tr>
      <tr>
        <td id="L437" class="blob-num js-line-number" data-line-number="437"></td>
        <td id="LC437" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">ml</span>[[<span class="pl-smi">k</span>]] <span class="pl-k">=</span> <span class="pl-smi">mConsist</span> <span class="pl-c"><span class="pl-c">#</span>initialize</span></td>
      </tr>
      <tr>
        <td id="L438" class="blob-num js-line-number" data-line-number="438"></td>
        <td id="LC438" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L439" class="blob-num js-line-number" data-line-number="439"></td>
        <td id="LC439" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">this_assignment</span><span class="pl-k">=</span><span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L440" class="blob-num js-line-number" data-line-number="440"></td>
        <td id="LC440" class="blob-code blob-code-inner js-file-line">      <span class="pl-k">if</span>(<span class="pl-smi">clusterAlg</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>hc<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L441" class="blob-num js-line-number" data-line-number="441"></td>
        <td id="LC441" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span>#prune to k for hc</span></td>
      </tr>
      <tr>
        <td id="L442" class="blob-num js-line-number" data-line-number="442"></td>
        <td id="LC442" class="blob-code blob-code-inner js-file-line">        <span class="pl-v">this_assignment</span> <span class="pl-k">=</span> cutree(<span class="pl-smi">this_cluster</span>,<span class="pl-smi">k</span>)</td>
      </tr>
      <tr>
        <td id="L443" class="blob-num js-line-number" data-line-number="443"></td>
        <td id="LC443" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L444" class="blob-num js-line-number" data-line-number="444"></td>
        <td id="LC444" class="blob-code blob-code-inner js-file-line">      }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">clusterAlg</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>kmdist<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L445" class="blob-num js-line-number" data-line-number="445"></td>
        <td id="LC445" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">this_assignment</span> <span class="pl-k">=</span> kmeans(<span class="pl-smi">this_dist</span>, <span class="pl-smi">k</span>, <span class="pl-v">iter.max</span> <span class="pl-k">=</span> <span class="pl-c1">10</span>, <span class="pl-v">nstart</span> <span class="pl-k">=</span> <span class="pl-c1">1</span>, <span class="pl-v">algorithm</span> <span class="pl-k">=</span> c(<span class="pl-s"><span class="pl-pds">&quot;</span>Hartigan-Wong<span class="pl-pds">&quot;</span></span>) )<span class="pl-k">$</span><span class="pl-smi">cluster</span></td>
      </tr>
      <tr>
        <td id="L446" class="blob-num js-line-number" data-line-number="446"></td>
        <td id="LC446" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L447" class="blob-num js-line-number" data-line-number="447"></td>
        <td id="LC447" class="blob-code blob-code-inner js-file-line">      }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">clusterAlg</span><span class="pl-k">==</span><span class="pl-s"><span class="pl-pds">&quot;</span>km<span class="pl-pds">&quot;</span></span>){</td>
      </tr>
      <tr>
        <td id="L448" class="blob-num js-line-number" data-line-number="448"></td>
        <td id="LC448" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span>#this_dist should now be a matrix corresponding to the result from sampleCols</span></td>
      </tr>
      <tr>
        <td id="L449" class="blob-num js-line-number" data-line-number="449"></td>
        <td id="LC449" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">this_assignment</span> <span class="pl-k">&lt;-</span> kmeans( t( <span class="pl-smi">this_dist</span> ),</td>
      </tr>
      <tr>
        <td id="L450" class="blob-num js-line-number" data-line-number="450"></td>
        <td id="LC450" class="blob-code blob-code-inner js-file-line">                                   <span class="pl-smi">k</span>,</td>
      </tr>
      <tr>
        <td id="L451" class="blob-num js-line-number" data-line-number="451"></td>
        <td id="LC451" class="blob-code blob-code-inner js-file-line">                                   <span class="pl-v">iter.max</span> <span class="pl-k">=</span> <span class="pl-c1">10</span>,</td>
      </tr>
      <tr>
        <td id="L452" class="blob-num js-line-number" data-line-number="452"></td>
        <td id="LC452" class="blob-code blob-code-inner js-file-line">                                   <span class="pl-v">nstart</span> <span class="pl-k">=</span> <span class="pl-c1">1</span>,</td>
      </tr>
      <tr>
        <td id="L453" class="blob-num js-line-number" data-line-number="453"></td>
        <td id="LC453" class="blob-code blob-code-inner js-file-line">                                   <span class="pl-v">algorithm</span> <span class="pl-k">=</span> c(<span class="pl-s"><span class="pl-pds">&quot;</span>Hartigan-Wong<span class="pl-pds">&quot;</span></span>) )<span class="pl-k">$</span><span class="pl-smi">cluster</span></td>
      </tr>
      <tr>
        <td id="L454" class="blob-num js-line-number" data-line-number="454"></td>
        <td id="LC454" class="blob-code blob-code-inner js-file-line">      }<span class="pl-k">else</span> <span class="pl-k">if</span> ( <span class="pl-smi">clusterAlg</span> <span class="pl-k">==</span> <span class="pl-s"><span class="pl-pds">&quot;</span>pam<span class="pl-pds">&quot;</span></span> ) {</td>
      </tr>
      <tr>
        <td id="L455" class="blob-num js-line-number" data-line-number="455"></td>
        <td id="LC455" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">this_assignment</span> <span class="pl-k">&lt;-</span> pam( <span class="pl-v">x</span><span class="pl-k">=</span><span class="pl-smi">this_dist</span>,</td>
      </tr>
      <tr>
        <td id="L456" class="blob-num js-line-number" data-line-number="456"></td>
        <td id="LC456" class="blob-code blob-code-inner js-file-line">                                <span class="pl-smi">k</span>,</td>
      </tr>
      <tr>
        <td id="L457" class="blob-num js-line-number" data-line-number="457"></td>
        <td id="LC457" class="blob-code blob-code-inner js-file-line">                                <span class="pl-v">diss</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span>,</td>
      </tr>
      <tr>
        <td id="L458" class="blob-num js-line-number" data-line-number="458"></td>
        <td id="LC458" class="blob-code blob-code-inner js-file-line">                                <span class="pl-v">metric</span><span class="pl-k">=</span><span class="pl-smi">distance</span>, </td>
      </tr>
      <tr>
        <td id="L459" class="blob-num js-line-number" data-line-number="459"></td>
        <td id="LC459" class="blob-code blob-code-inner js-file-line">                                <span class="pl-v">cluster.only</span><span class="pl-k">=</span><span class="pl-c1">TRUE</span> )</td>
      </tr>
      <tr>
        <td id="L460" class="blob-num js-line-number" data-line-number="460"></td>
        <td id="LC460" class="blob-code blob-code-inner js-file-line">      } <span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L461" class="blob-num js-line-number" data-line-number="461"></td>
        <td id="LC461" class="blob-code blob-code-inner js-file-line">        <span class="pl-c"><span class="pl-c">#</span>#optional cluterArg Hook.</span></td>
      </tr>
      <tr>
        <td id="L462" class="blob-num js-line-number" data-line-number="462"></td>
        <td id="LC462" class="blob-code blob-code-inner js-file-line">        <span class="pl-smi">this_assignment</span> <span class="pl-k">&lt;-</span> get(<span class="pl-smi">clusterAlg</span>)(<span class="pl-smi">this_dist</span>, <span class="pl-smi">k</span>)</td>
      </tr>
      <tr>
        <td id="L463" class="blob-num js-line-number" data-line-number="463"></td>
        <td id="LC463" class="blob-code blob-code-inner js-file-line">      }</td>
      </tr>
      <tr>
        <td id="L464" class="blob-num js-line-number" data-line-number="464"></td>
        <td id="LC464" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span>#add to tally				</span></td>
      </tr>
      <tr>
        <td id="L465" class="blob-num js-line-number" data-line-number="465"></td>
        <td id="LC465" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">ml</span>[[<span class="pl-smi">k</span>]] <span class="pl-k">&lt;-</span> connectivityMatrix( <span class="pl-smi">this_assignment</span>,</td>
      </tr>
      <tr>
        <td id="L466" class="blob-num js-line-number" data-line-number="466"></td>
        <td id="LC466" class="blob-code blob-code-inner js-file-line">                                     <span class="pl-smi">ml</span>[[<span class="pl-smi">k</span>]],</td>
      </tr>
      <tr>
        <td id="L467" class="blob-num js-line-number" data-line-number="467"></td>
        <td id="LC467" class="blob-code blob-code-inner js-file-line">                                     <span class="pl-smi">sample_x</span>[[<span class="pl-c1">3</span>]] )</td>
      </tr>
      <tr>
        <td id="L468" class="blob-num js-line-number" data-line-number="468"></td>
        <td id="LC468" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L469" class="blob-num js-line-number" data-line-number="469"></td>
        <td id="LC469" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L470" class="blob-num js-line-number" data-line-number="470"></td>
        <td id="LC470" class="blob-code blob-code-inner js-file-line">	</td>
      </tr>
      <tr>
        <td id="L471" class="blob-num js-line-number" data-line-number="471"></td>
        <td id="LC471" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L472" class="blob-num js-line-number" data-line-number="472"></td>
        <td id="LC472" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#consensus fraction</span></td>
      </tr>
      <tr>
        <td id="L473" class="blob-num js-line-number" data-line-number="473"></td>
        <td id="LC473" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">res</span> <span class="pl-k">=</span> vector(<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>list<span class="pl-pds">&quot;</span></span>,<span class="pl-smi">maxK</span>)</td>
      </tr>
      <tr>
        <td id="L474" class="blob-num js-line-number" data-line-number="474"></td>
        <td id="LC474" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> (<span class="pl-smi">k</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span><span class="pl-smi">maxK</span>){</td>
      </tr>
      <tr>
        <td id="L475" class="blob-num js-line-number" data-line-number="475"></td>
        <td id="LC475" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>#fill in other half of matrix for tally and count.</span></td>
      </tr>
      <tr>
        <td id="L476" class="blob-num js-line-number" data-line-number="476"></td>
        <td id="LC476" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">tmp</span> <span class="pl-k">=</span> triangle(<span class="pl-smi">ml</span>[[<span class="pl-smi">k</span>]],<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-c1">3</span>)</td>
      </tr>
      <tr>
        <td id="L477" class="blob-num js-line-number" data-line-number="477"></td>
        <td id="LC477" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">tmpCount</span> <span class="pl-k">=</span> triangle(<span class="pl-smi">mCount</span>,<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-c1">3</span>)</td>
      </tr>
      <tr>
        <td id="L478" class="blob-num js-line-number" data-line-number="478"></td>
        <td id="LC478" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">res</span>[[<span class="pl-smi">k</span>]] <span class="pl-k">=</span> <span class="pl-smi">tmp</span> <span class="pl-k">/</span> <span class="pl-smi">tmpCount</span></td>
      </tr>
      <tr>
        <td id="L479" class="blob-num js-line-number" data-line-number="479"></td>
        <td id="LC479" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">res</span>[[<span class="pl-smi">k</span>]][which(<span class="pl-smi">tmpCount</span><span class="pl-k">==</span><span class="pl-c1">0</span>)] <span class="pl-k">=</span> <span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L480" class="blob-num js-line-number" data-line-number="480"></td>
        <td id="LC480" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L481" class="blob-num js-line-number" data-line-number="481"></td>
        <td id="LC481" class="blob-code blob-code-inner js-file-line">  message(<span class="pl-s"><span class="pl-pds">&quot;</span>end fraction<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L482" class="blob-num js-line-number" data-line-number="482"></td>
        <td id="LC482" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">return</span>(<span class="pl-smi">res</span>)</td>
      </tr>
      <tr>
        <td id="L483" class="blob-num js-line-number" data-line-number="483"></td>
        <td id="LC483" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L484" class="blob-num js-line-number" data-line-number="484"></td>
        <td id="LC484" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L485" class="blob-num js-line-number" data-line-number="485"></td>
        <td id="LC485" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L486" class="blob-num js-line-number" data-line-number="486"></td>
        <td id="LC486" class="blob-code blob-code-inner js-file-line"><span class="pl-en">connectivityMatrix</span> <span class="pl-k">&lt;-</span> <span class="pl-k">function</span>( <span class="pl-smi">clusterAssignments</span>, <span class="pl-smi">m</span>, <span class="pl-smi">sampleKey</span>){</td>
      </tr>
      <tr>
        <td id="L487" class="blob-num js-line-number" data-line-number="487"></td>
        <td id="LC487" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#input: named vector of cluster assignments, matrix to add connectivities</span></td>
      </tr>
      <tr>
        <td id="L488" class="blob-num js-line-number" data-line-number="488"></td>
        <td id="LC488" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>#output: connectivity matrix</span></td>
      </tr>
      <tr>
        <td id="L489" class="blob-num js-line-number" data-line-number="489"></td>
        <td id="LC489" class="blob-code blob-code-inner js-file-line">  names( <span class="pl-smi">clusterAssignments</span> ) <span class="pl-k">&lt;-</span> <span class="pl-smi">sampleKey</span> </td>
      </tr>
      <tr>
        <td id="L490" class="blob-num js-line-number" data-line-number="490"></td>
        <td id="LC490" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">cls</span> <span class="pl-k">&lt;-</span> lapply( unique( <span class="pl-smi">clusterAssignments</span> ), <span class="pl-k">function</span>(<span class="pl-smi">i</span>) as.numeric( names( <span class="pl-smi">clusterAssignments</span>[ <span class="pl-smi">clusterAssignments</span> <span class="pl-k">%in%</span> <span class="pl-smi">i</span> ] ) ) )  <span class="pl-c"><span class="pl-c">#</span>list samples by clusterId</span></td>
      </tr>
      <tr>
        <td id="L491" class="blob-num js-line-number" data-line-number="491"></td>
        <td id="LC491" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L492" class="blob-num js-line-number" data-line-number="492"></td>
        <td id="LC492" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> ( <span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span>length( <span class="pl-smi">cls</span> ) ) {</td>
      </tr>
      <tr>
        <td id="L493" class="blob-num js-line-number" data-line-number="493"></td>
        <td id="LC493" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">nelts</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">1</span><span class="pl-k">:</span>ncol( <span class="pl-smi">m</span> )</td>
      </tr>
      <tr>
        <td id="L494" class="blob-num js-line-number" data-line-number="494"></td>
        <td id="LC494" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">cl</span> <span class="pl-k">&lt;-</span> as.numeric( <span class="pl-smi">nelts</span> <span class="pl-k">%in%</span> <span class="pl-smi">cls</span>[[<span class="pl-smi">i</span>]] ) <span class="pl-c"><span class="pl-c">#</span># produces a binary vector</span></td>
      </tr>
      <tr>
        <td id="L495" class="blob-num js-line-number" data-line-number="495"></td>
        <td id="LC495" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">updt</span> <span class="pl-k">&lt;-</span> outer( <span class="pl-smi">cl</span>, <span class="pl-smi">cl</span> ) <span class="pl-c"><span class="pl-c">#</span>product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;</span></td>
      </tr>
      <tr>
        <td id="L496" class="blob-num js-line-number" data-line-number="496"></td>
        <td id="LC496" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">m</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">m</span> <span class="pl-k">+</span> <span class="pl-smi">updt</span></td>
      </tr>
      <tr>
        <td id="L497" class="blob-num js-line-number" data-line-number="497"></td>
        <td id="LC497" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L498" class="blob-num js-line-number" data-line-number="498"></td>
        <td id="LC498" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">return</span>(<span class="pl-smi">m</span>)</td>
      </tr>
      <tr>
        <td id="L499" class="blob-num js-line-number" data-line-number="499"></td>
        <td id="LC499" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L500" class="blob-num js-line-number" data-line-number="500"></td>
        <td id="LC500" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L501" class="blob-num js-line-number" data-line-number="501"></td>
        <td id="LC501" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L502" class="blob-num js-line-number" data-line-number="502"></td>
        <td id="LC502" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L503" class="blob-num js-line-number" data-line-number="503"></td>
        <td id="LC503" class="blob-code blob-code-inner js-file-line"><span class="pl-en">sampleCols</span> <span class="pl-k">&lt;-</span> <span class="pl-k">function</span>( <span class="pl-smi">d</span>,</td>
      </tr>
      <tr>
        <td id="L504" class="blob-num js-line-number" data-line-number="504"></td>
        <td id="LC504" class="blob-code blob-code-inner js-file-line">                        <span class="pl-v">pSamp</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L505" class="blob-num js-line-number" data-line-number="505"></td>
        <td id="LC505" class="blob-code blob-code-inner js-file-line">                        <span class="pl-v">pRow</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L506" class="blob-num js-line-number" data-line-number="506"></td>
        <td id="LC506" class="blob-code blob-code-inner js-file-line">                        <span class="pl-v">weightsItem</span><span class="pl-k">=</span><span class="pl-c1">NULL</span>,</td>
      </tr>
      <tr>
        <td id="L507" class="blob-num js-line-number" data-line-number="507"></td>
        <td id="LC507" class="blob-code blob-code-inner js-file-line">                        <span class="pl-v">weightsFeature</span><span class="pl-k">=</span><span class="pl-c1">NULL</span> ){</td>
      </tr>
      <tr>
        <td id="L508" class="blob-num js-line-number" data-line-number="508"></td>
        <td id="LC508" class="blob-code blob-code-inner js-file-line"><span class="pl-c"><span class="pl-c">#</span># returns a list with the sample columns, as well as the sub-matrix &amp; sample features (if necessary)</span></td>
      </tr>
      <tr>
        <td id="L509" class="blob-num js-line-number" data-line-number="509"></td>
        <td id="LC509" class="blob-code blob-code-inner js-file-line"><span class="pl-c"><span class="pl-c">#</span>#  if no sampling over the features is performed, the submatrix &amp; sample features are returned as NAs</span></td>
      </tr>
      <tr>
        <td id="L510" class="blob-num js-line-number" data-line-number="510"></td>
        <td id="LC510" class="blob-code blob-code-inner js-file-line"><span class="pl-c"><span class="pl-c">#</span>#  to reduce memory overhead</span></td>
      </tr>
      <tr>
        <td id="L511" class="blob-num js-line-number" data-line-number="511"></td>
        <td id="LC511" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L512" class="blob-num js-line-number" data-line-number="512"></td>
        <td id="LC512" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L513" class="blob-num js-line-number" data-line-number="513"></td>
        <td id="LC513" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">space</span> <span class="pl-k">&lt;-</span> ifelse( inherits( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>dist<span class="pl-pds">&quot;</span></span> ), ncol( as.matrix(<span class="pl-smi">d</span>) ), ncol(<span class="pl-smi">d</span>) )</td>
      </tr>
      <tr>
        <td id="L514" class="blob-num js-line-number" data-line-number="514"></td>
        <td id="LC514" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">sampleN</span> <span class="pl-k">&lt;-</span> floor(<span class="pl-smi">space</span><span class="pl-k">*</span><span class="pl-smi">pSamp</span>)</td>
      </tr>
      <tr>
        <td id="L515" class="blob-num js-line-number" data-line-number="515"></td>
        <td id="LC515" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">sampCols</span> <span class="pl-k">&lt;-</span> sort( sample(<span class="pl-smi">space</span>, <span class="pl-smi">sampleN</span>, <span class="pl-v">replace</span> <span class="pl-k">=</span> <span class="pl-c1">FALSE</span>, <span class="pl-v">prob</span> <span class="pl-k">=</span> <span class="pl-smi">weightsItem</span>) )</td>
      </tr>
      <tr>
        <td id="L516" class="blob-num js-line-number" data-line-number="516"></td>
        <td id="LC516" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L517" class="blob-num js-line-number" data-line-number="517"></td>
        <td id="LC517" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">this_sample</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">sampRows</span> <span class="pl-k">&lt;-</span> <span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L518" class="blob-num js-line-number" data-line-number="518"></td>
        <td id="LC518" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span> ( inherits( <span class="pl-smi">d</span>, <span class="pl-s"><span class="pl-pds">&quot;</span>matrix<span class="pl-pds">&quot;</span></span> ) ) {</td>
      </tr>
      <tr>
        <td id="L519" class="blob-num js-line-number" data-line-number="519"></td>
        <td id="LC519" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">if</span> ( (<span class="pl-k">!</span> is.null( <span class="pl-smi">pRow</span> ) ) <span class="pl-k">&amp;&amp;</span></td>
      </tr>
      <tr>
        <td id="L520" class="blob-num js-line-number" data-line-number="520"></td>
        <td id="LC520" class="blob-code blob-code-inner js-file-line">         ( (<span class="pl-smi">pRow</span> <span class="pl-k">&lt;</span> <span class="pl-c1">1</span> ) <span class="pl-k">||</span> (<span class="pl-k">!</span> is.null( <span class="pl-smi">weightsFeature</span> ) ) ) ) {</td>
      </tr>
      <tr>
        <td id="L521" class="blob-num js-line-number" data-line-number="521"></td>
        <td id="LC521" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># only sample the rows and generate a sub-matrix if we&#39;re sampling over the row/gene/features</span></td>
      </tr>
      <tr>
        <td id="L522" class="blob-num js-line-number" data-line-number="522"></td>
        <td id="LC522" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">space</span> <span class="pl-k">=</span> nrow(<span class="pl-smi">d</span>)</td>
      </tr>
      <tr>
        <td id="L523" class="blob-num js-line-number" data-line-number="523"></td>
        <td id="LC523" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">sampleN</span> <span class="pl-k">=</span> floor(<span class="pl-smi">space</span><span class="pl-k">*</span><span class="pl-smi">pRow</span>)</td>
      </tr>
      <tr>
        <td id="L524" class="blob-num js-line-number" data-line-number="524"></td>
        <td id="LC524" class="blob-code blob-code-inner js-file-line">      <span class="pl-v">sampRows</span> <span class="pl-k">=</span> sort( sample(<span class="pl-smi">space</span>, <span class="pl-smi">sampleN</span>, <span class="pl-v">replace</span> <span class="pl-k">=</span> <span class="pl-c1">FALSE</span>, <span class="pl-v">prob</span> <span class="pl-k">=</span> <span class="pl-smi">weightsFeature</span>) )</td>
      </tr>
      <tr>
        <td id="L525" class="blob-num js-line-number" data-line-number="525"></td>
        <td id="LC525" class="blob-code blob-code-inner js-file-line">      <span class="pl-smi">this_sample</span> <span class="pl-k">&lt;-</span> <span class="pl-smi">d</span>[<span class="pl-smi">sampRows</span>,<span class="pl-smi">sampCols</span>]</td>
      </tr>
      <tr>
        <td id="L526" class="blob-num js-line-number" data-line-number="526"></td>
        <td id="LC526" class="blob-code blob-code-inner js-file-line">      dimnames(<span class="pl-smi">this_sample</span>) <span class="pl-k">&lt;-</span> <span class="pl-c1">NULL</span></td>
      </tr>
      <tr>
        <td id="L527" class="blob-num js-line-number" data-line-number="527"></td>
        <td id="LC527" class="blob-code blob-code-inner js-file-line">    } <span class="pl-k">else</span> {</td>
      </tr>
      <tr>
        <td id="L528" class="blob-num js-line-number" data-line-number="528"></td>
        <td id="LC528" class="blob-code blob-code-inner js-file-line">      <span class="pl-c"><span class="pl-c">#</span># do nothing</span></td>
      </tr>
      <tr>
        <td id="L529" class="blob-num js-line-number" data-line-number="529"></td>
        <td id="LC529" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L530" class="blob-num js-line-number" data-line-number="530"></td>
        <td id="LC530" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L531" class="blob-num js-line-number" data-line-number="531"></td>
        <td id="LC531" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">return</span>( <span class="pl-k">list</span>( <span class="pl-v">submat</span><span class="pl-k">=</span><span class="pl-smi">this_sample</span>,</td>
      </tr>
      <tr>
        <td id="L532" class="blob-num js-line-number" data-line-number="532"></td>
        <td id="LC532" class="blob-code blob-code-inner js-file-line">                <span class="pl-v">subrows</span><span class="pl-k">=</span><span class="pl-smi">sampRows</span>,</td>
      </tr>
      <tr>
        <td id="L533" class="blob-num js-line-number" data-line-number="533"></td>
        <td id="LC533" class="blob-code blob-code-inner js-file-line">                <span class="pl-v">subcols</span><span class="pl-k">=</span><span class="pl-smi">sampCols</span> ) )</td>
      </tr>
      <tr>
        <td id="L534" class="blob-num js-line-number" data-line-number="534"></td>
        <td id="LC534" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L535" class="blob-num js-line-number" data-line-number="535"></td>
        <td id="LC535" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L536" class="blob-num js-line-number" data-line-number="536"></td>
        <td id="LC536" class="blob-code blob-code-inner js-file-line"><span class="pl-v">CDF</span><span class="pl-k">=</span><span class="pl-k">function</span>(<span class="pl-smi">ml</span>,<span class="pl-v">breaks</span><span class="pl-k">=</span><span class="pl-c1">100</span>){</td>
      </tr>
      <tr>
        <td id="L537" class="blob-num js-line-number" data-line-number="537"></td>
        <td id="LC537" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>plot CDF distribution</span></td>
      </tr>
      <tr>
        <td id="L538" class="blob-num js-line-number" data-line-number="538"></td>
        <td id="LC538" class="blob-code blob-code-inner js-file-line">  plot(c(<span class="pl-c1">0</span>),<span class="pl-v">xlim</span><span class="pl-k">=</span>c(<span class="pl-c1">0</span>,<span class="pl-c1">1</span>),<span class="pl-v">ylim</span><span class="pl-k">=</span>c(<span class="pl-c1">0</span>,<span class="pl-c1">1</span>),<span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span>,<span class="pl-v">bg</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>white<span class="pl-pds">&quot;</span></span>,<span class="pl-v">xlab</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>consensus index<span class="pl-pds">&quot;</span></span>,<span class="pl-v">ylab</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>CDF<span class="pl-pds">&quot;</span></span>,<span class="pl-v">main</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>consensus CDF<span class="pl-pds">&quot;</span></span>, <span class="pl-v">las</span><span class="pl-k">=</span><span class="pl-c1">2</span>)</td>
      </tr>
      <tr>
        <td id="L539" class="blob-num js-line-number" data-line-number="539"></td>
        <td id="LC539" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">k</span><span class="pl-k">=</span>length(<span class="pl-smi">ml</span>)</td>
      </tr>
      <tr>
        <td id="L540" class="blob-num js-line-number" data-line-number="540"></td>
        <td id="LC540" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">this_colors</span> <span class="pl-k">=</span> rainbow(<span class="pl-smi">k</span><span class="pl-k">-</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L541" class="blob-num js-line-number" data-line-number="541"></td>
        <td id="LC541" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">areaK</span> <span class="pl-k">=</span> c()</td>
      </tr>
      <tr>
        <td id="L542" class="blob-num js-line-number" data-line-number="542"></td>
        <td id="LC542" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span> (<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span>length(<span class="pl-smi">ml</span>)){</td>
      </tr>
      <tr>
        <td id="L543" class="blob-num js-line-number" data-line-number="543"></td>
        <td id="LC543" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">v</span><span class="pl-k">=</span>triangle(<span class="pl-smi">ml</span>[[<span class="pl-smi">i</span>]],<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-c1">1</span>)</td>
      </tr>
      <tr>
        <td id="L544" class="blob-num js-line-number" data-line-number="544"></td>
        <td id="LC544" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L545" class="blob-num js-line-number" data-line-number="545"></td>
        <td id="LC545" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>empirical CDF distribution. default number of breaks is 100    </span></td>
      </tr>
      <tr>
        <td id="L546" class="blob-num js-line-number" data-line-number="546"></td>
        <td id="LC546" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">h</span> <span class="pl-k">=</span> hist(<span class="pl-smi">v</span>, <span class="pl-v">plot</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span>, <span class="pl-v">breaks</span><span class="pl-k">=</span>seq(<span class="pl-c1">0</span>,<span class="pl-c1">1</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span><span class="pl-k">/</span><span class="pl-smi">breaks</span>))</td>
      </tr>
      <tr>
        <td id="L547" class="blob-num js-line-number" data-line-number="547"></td>
        <td id="LC547" class="blob-code blob-code-inner js-file-line">    <span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-v">counts</span> <span class="pl-k">=</span> cumsum(<span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">counts</span>)<span class="pl-k">/</span>sum(<span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">counts</span>)</td>
      </tr>
      <tr>
        <td id="L548" class="blob-num js-line-number" data-line-number="548"></td>
        <td id="LC548" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L549" class="blob-num js-line-number" data-line-number="549"></td>
        <td id="LC549" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>calculate area under CDF curve, by histogram method.</span></td>
      </tr>
      <tr>
        <td id="L550" class="blob-num js-line-number" data-line-number="550"></td>
        <td id="LC550" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">thisArea</span><span class="pl-k">=</span><span class="pl-c1">0</span></td>
      </tr>
      <tr>
        <td id="L551" class="blob-num js-line-number" data-line-number="551"></td>
        <td id="LC551" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">for</span> (<span class="pl-smi">bi</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span>(length(<span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">breaks</span>)<span class="pl-k">-</span><span class="pl-c1">1</span>)){</td>
      </tr>
      <tr>
        <td id="L552" class="blob-num js-line-number" data-line-number="552"></td>
        <td id="LC552" class="blob-code blob-code-inner js-file-line">       <span class="pl-v">thisArea</span> <span class="pl-k">=</span> <span class="pl-smi">thisArea</span> <span class="pl-k">+</span> <span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">counts</span>[<span class="pl-smi">bi</span>]<span class="pl-k">*</span>(<span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">breaks</span>[<span class="pl-smi">bi</span><span class="pl-k">+</span><span class="pl-c1">1</span>]<span class="pl-k">-</span><span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">breaks</span>[<span class="pl-smi">bi</span>]) <span class="pl-c"><span class="pl-c">#</span>increment by height by width</span></td>
      </tr>
      <tr>
        <td id="L553" class="blob-num js-line-number" data-line-number="553"></td>
        <td id="LC553" class="blob-code blob-code-inner js-file-line">       <span class="pl-v">bi</span> <span class="pl-k">=</span> <span class="pl-smi">bi</span> <span class="pl-k">+</span> <span class="pl-c1">1</span></td>
      </tr>
      <tr>
        <td id="L554" class="blob-num js-line-number" data-line-number="554"></td>
        <td id="LC554" class="blob-code blob-code-inner js-file-line">    }</td>
      </tr>
      <tr>
        <td id="L555" class="blob-num js-line-number" data-line-number="555"></td>
        <td id="LC555" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">areaK</span> <span class="pl-k">=</span> c(<span class="pl-smi">areaK</span>,<span class="pl-smi">thisArea</span>)</td>
      </tr>
      <tr>
        <td id="L556" class="blob-num js-line-number" data-line-number="556"></td>
        <td id="LC556" class="blob-code blob-code-inner js-file-line">    lines(<span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">mids</span>,<span class="pl-smi">h</span><span class="pl-k">$</span><span class="pl-smi">counts</span>,<span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">this_colors</span>[<span class="pl-smi">i</span><span class="pl-k">-</span><span class="pl-c1">1</span>],<span class="pl-v">lwd</span><span class="pl-k">=</span><span class="pl-c1">2</span>,<span class="pl-v">type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&#39;</span>l<span class="pl-pds">&#39;</span></span>)</td>
      </tr>
      <tr>
        <td id="L557" class="blob-num js-line-number" data-line-number="557"></td>
        <td id="LC557" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L558" class="blob-num js-line-number" data-line-number="558"></td>
        <td id="LC558" class="blob-code blob-code-inner js-file-line">  legend(<span class="pl-c1">0.8</span>,<span class="pl-c1">0.5</span>,<span class="pl-v">legend</span><span class="pl-k">=</span>paste(rep(<span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>,<span class="pl-smi">k</span><span class="pl-k">-</span><span class="pl-c1">1</span>),seq(<span class="pl-c1">2</span>,<span class="pl-smi">k</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span>),<span class="pl-v">sep</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span><span class="pl-pds">&quot;</span></span>),<span class="pl-v">fill</span><span class="pl-k">=</span><span class="pl-smi">this_colors</span>)</td>
      </tr>
      <tr>
        <td id="L559" class="blob-num js-line-number" data-line-number="559"></td>
        <td id="LC559" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L560" class="blob-num js-line-number" data-line-number="560"></td>
        <td id="LC560" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>plot area under CDF change.</span></td>
      </tr>
      <tr>
        <td id="L561" class="blob-num js-line-number" data-line-number="561"></td>
        <td id="LC561" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">deltaK</span><span class="pl-k">=</span><span class="pl-smi">areaK</span>[<span class="pl-c1">1</span>] <span class="pl-c"><span class="pl-c">#</span>initial auc at k=2</span></td>
      </tr>
      <tr>
        <td id="L562" class="blob-num js-line-number" data-line-number="562"></td>
        <td id="LC562" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span>(<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span>(length(<span class="pl-smi">areaK</span>))){</td>
      </tr>
      <tr>
        <td id="L563" class="blob-num js-line-number" data-line-number="563"></td>
        <td id="LC563" class="blob-code blob-code-inner js-file-line">    <span class="pl-c"><span class="pl-c">#</span>proportional increase relative to prior K.</span></td>
      </tr>
      <tr>
        <td id="L564" class="blob-num js-line-number" data-line-number="564"></td>
        <td id="LC564" class="blob-code blob-code-inner js-file-line">    <span class="pl-v">deltaK</span> <span class="pl-k">=</span> c(<span class="pl-smi">deltaK</span>,( <span class="pl-smi">areaK</span>[<span class="pl-smi">i</span>] <span class="pl-k">-</span> <span class="pl-smi">areaK</span>[<span class="pl-smi">i</span><span class="pl-k">-</span><span class="pl-c1">1</span>])<span class="pl-k">/</span><span class="pl-smi">areaK</span>[<span class="pl-smi">i</span><span class="pl-k">-</span><span class="pl-c1">1</span>])</td>
      </tr>
      <tr>
        <td id="L565" class="blob-num js-line-number" data-line-number="565"></td>
        <td id="LC565" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L566" class="blob-num js-line-number" data-line-number="566"></td>
        <td id="LC566" class="blob-code blob-code-inner js-file-line">  plot(<span class="pl-c1">1</span><span class="pl-k">+</span>(<span class="pl-c1">1</span><span class="pl-k">:</span>length(<span class="pl-smi">deltaK</span>)),<span class="pl-v">y</span><span class="pl-k">=</span><span class="pl-smi">deltaK</span>,<span class="pl-v">xlab</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>k<span class="pl-pds">&quot;</span></span>,<span class="pl-v">ylab</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>relative change in area under CDF curve<span class="pl-pds">&quot;</span></span>,<span class="pl-v">main</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>Delta area<span class="pl-pds">&quot;</span></span>,<span class="pl-v">type</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>b<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L567" class="blob-num js-line-number" data-line-number="567"></td>
        <td id="LC567" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L568" class="blob-num js-line-number" data-line-number="568"></td>
        <td id="LC568" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L569" class="blob-num js-line-number" data-line-number="569"></td>
        <td id="LC569" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L570" class="blob-num js-line-number" data-line-number="570"></td>
        <td id="LC570" class="blob-code blob-code-inner js-file-line"><span class="pl-v">myPal</span> <span class="pl-k">=</span> <span class="pl-k">function</span>(<span class="pl-v">n</span><span class="pl-k">=</span><span class="pl-c1">10</span>){</td>
      </tr>
      <tr>
        <td id="L571" class="blob-num js-line-number" data-line-number="571"></td>
        <td id="LC571" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>returns n colors</span></td>
      </tr>
      <tr>
        <td id="L572" class="blob-num js-line-number" data-line-number="572"></td>
        <td id="LC572" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">seq</span> <span class="pl-k">=</span> rev(seq(<span class="pl-c1">0</span>,<span class="pl-c1">255</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">255</span><span class="pl-k">/</span>(<span class="pl-smi">n</span>)))</td>
      </tr>
      <tr>
        <td id="L573" class="blob-num js-line-number" data-line-number="573"></td>
        <td id="LC573" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">palRGB</span> <span class="pl-k">=</span> cbind(<span class="pl-smi">seq</span>,<span class="pl-smi">seq</span>,<span class="pl-c1">255</span>)</td>
      </tr>
      <tr>
        <td id="L574" class="blob-num js-line-number" data-line-number="574"></td>
        <td id="LC574" class="blob-code blob-code-inner js-file-line">  rgb(<span class="pl-smi">palRGB</span>,<span class="pl-v">maxColorValue</span><span class="pl-k">=</span><span class="pl-c1">255</span>)</td>
      </tr>
      <tr>
        <td id="L575" class="blob-num js-line-number" data-line-number="575"></td>
        <td id="LC575" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L576" class="blob-num js-line-number" data-line-number="576"></td>
        <td id="LC576" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L577" class="blob-num js-line-number" data-line-number="577"></td>
        <td id="LC577" class="blob-code blob-code-inner js-file-line"><span class="pl-v">setClusterColors</span> <span class="pl-k">=</span> <span class="pl-k">function</span>(<span class="pl-smi">past_ct</span>,<span class="pl-smi">ct</span>,<span class="pl-smi">colorU</span>,<span class="pl-smi">colorList</span>){</td>
      </tr>
      <tr>
        <td id="L578" class="blob-num js-line-number" data-line-number="578"></td>
        <td id="LC578" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">#</span>description: sets common color of clusters between different K</span></td>
      </tr>
      <tr>
        <td id="L579" class="blob-num js-line-number" data-line-number="579"></td>
        <td id="LC579" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">newColors</span> <span class="pl-k">=</span> c()</td>
      </tr>
      <tr>
        <td id="L580" class="blob-num js-line-number" data-line-number="580"></td>
        <td id="LC580" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">if</span>(length(<span class="pl-smi">colorList</span>)<span class="pl-k">==</span><span class="pl-c1">0</span>){</td>
      </tr>
      <tr>
        <td id="L581" class="blob-num js-line-number" data-line-number="581"></td>
        <td id="LC581" class="blob-code blob-code-inner js-file-line">		<span class="pl-c"><span class="pl-c">#</span>k==2</span></td>
      </tr>
      <tr>
        <td id="L582" class="blob-num js-line-number" data-line-number="582"></td>
        <td id="LC582" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">newColors</span> <span class="pl-k">=</span> <span class="pl-smi">colorU</span>[<span class="pl-smi">ct</span>]</td>
      </tr>
      <tr>
        <td id="L583" class="blob-num js-line-number" data-line-number="583"></td>
        <td id="LC583" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">colori</span><span class="pl-k">=</span><span class="pl-c1">2</span></td>
      </tr>
      <tr>
        <td id="L584" class="blob-num js-line-number" data-line-number="584"></td>
        <td id="LC584" class="blob-code blob-code-inner js-file-line">	}<span class="pl-k">else</span>{</td>
      </tr>
      <tr>
        <td id="L585" class="blob-num js-line-number" data-line-number="585"></td>
        <td id="LC585" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">newColors</span> <span class="pl-k">=</span> rep(<span class="pl-c1">NULL</span>,length(<span class="pl-smi">ct</span>))</td>
      </tr>
      <tr>
        <td id="L586" class="blob-num js-line-number" data-line-number="586"></td>
        <td id="LC586" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">colori</span> <span class="pl-k">=</span> <span class="pl-smi">colorList</span>[[<span class="pl-c1">2</span>]]</td>
      </tr>
      <tr>
        <td id="L587" class="blob-num js-line-number" data-line-number="587"></td>
        <td id="LC587" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">mo</span><span class="pl-k">=</span>table(<span class="pl-smi">past_ct</span>,<span class="pl-smi">ct</span>)</td>
      </tr>
      <tr>
        <td id="L588" class="blob-num js-line-number" data-line-number="588"></td>
        <td id="LC588" class="blob-code blob-code-inner js-file-line">		<span class="pl-v">m</span><span class="pl-k">=</span><span class="pl-smi">mo</span><span class="pl-k">/</span>apply(<span class="pl-smi">mo</span>,<span class="pl-c1">1</span>,<span class="pl-smi">sum</span>)</td>
      </tr>
      <tr>
        <td id="L589" class="blob-num js-line-number" data-line-number="589"></td>
        <td id="LC589" class="blob-code blob-code-inner js-file-line">			<span class="pl-k">for</span>(<span class="pl-smi">tci</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span>ncol(<span class="pl-smi">m</span>)){ <span class="pl-c"><span class="pl-c">#</span> for each cluster</span></td>
      </tr>
      <tr>
        <td id="L590" class="blob-num js-line-number" data-line-number="590"></td>
        <td id="LC590" class="blob-code blob-code-inner js-file-line">				<span class="pl-v">maxC</span> <span class="pl-k">=</span> max(<span class="pl-smi">m</span>[,<span class="pl-smi">tci</span>])</td>
      </tr>
      <tr>
        <td id="L591" class="blob-num js-line-number" data-line-number="591"></td>
        <td id="LC591" class="blob-code blob-code-inner js-file-line">				<span class="pl-v">pci</span> <span class="pl-k">=</span> which(<span class="pl-smi">m</span>[,<span class="pl-smi">tci</span>] <span class="pl-k">==</span> <span class="pl-smi">maxC</span>)				</td>
      </tr>
      <tr>
        <td id="L592" class="blob-num js-line-number" data-line-number="592"></td>
        <td id="LC592" class="blob-code blob-code-inner js-file-line">				<span class="pl-k">if</span>( sum(<span class="pl-smi">m</span>[,<span class="pl-smi">tci</span>]<span class="pl-k">==</span><span class="pl-smi">maxC</span>)<span class="pl-k">==</span><span class="pl-c1">1</span> <span class="pl-k">&amp;</span> max(<span class="pl-smi">m</span>[<span class="pl-smi">pci</span>,])<span class="pl-k">==</span><span class="pl-smi">maxC</span> <span class="pl-k">&amp;</span> sum(<span class="pl-smi">m</span>[<span class="pl-smi">pci</span>,]<span class="pl-k">==</span><span class="pl-smi">maxC</span>)<span class="pl-k">==</span><span class="pl-c1">1</span>  )  {</td>
      </tr>
      <tr>
        <td id="L593" class="blob-num js-line-number" data-line-number="593"></td>
        <td id="LC593" class="blob-code blob-code-inner js-file-line">				<span class="pl-c"><span class="pl-c">#</span>if new column maximum is unique, same cell is row maximum and is also unique</span></td>
      </tr>
      <tr>
        <td id="L594" class="blob-num js-line-number" data-line-number="594"></td>
        <td id="LC594" class="blob-code blob-code-inner js-file-line">				<span class="pl-c"><span class="pl-c">#</span>#Note: the greatest of the prior clusters&#39; members are the greatest in a current cluster&#39;s members.</span></td>
      </tr>
      <tr>
        <td id="L595" class="blob-num js-line-number" data-line-number="595"></td>
        <td id="LC595" class="blob-code blob-code-inner js-file-line">					<span class="pl-smi">newColors</span>[which(<span class="pl-smi">ct</span><span class="pl-k">==</span><span class="pl-smi">tci</span>)] <span class="pl-k">=</span> unique(<span class="pl-smi">colorList</span>[[<span class="pl-c1">1</span>]][which(<span class="pl-smi">past_ct</span><span class="pl-k">==</span><span class="pl-smi">pci</span>)]) <span class="pl-c"><span class="pl-c">#</span> one value</span></td>
      </tr>
      <tr>
        <td id="L596" class="blob-num js-line-number" data-line-number="596"></td>
        <td id="LC596" class="blob-code blob-code-inner js-file-line">				}<span class="pl-k">else</span>{ <span class="pl-c"><span class="pl-c">#</span>add new color.</span></td>
      </tr>
      <tr>
        <td id="L597" class="blob-num js-line-number" data-line-number="597"></td>
        <td id="LC597" class="blob-code blob-code-inner js-file-line">					<span class="pl-v">colori</span><span class="pl-k">=</span><span class="pl-smi">colori</span><span class="pl-k">+</span><span class="pl-c1">1</span></td>
      </tr>
      <tr>
        <td id="L598" class="blob-num js-line-number" data-line-number="598"></td>
        <td id="LC598" class="blob-code blob-code-inner js-file-line">					<span class="pl-smi">newColors</span>[which(<span class="pl-smi">ct</span><span class="pl-k">==</span><span class="pl-smi">tci</span>)] <span class="pl-k">=</span> <span class="pl-smi">colorU</span>[<span class="pl-smi">colori</span>]</td>
      </tr>
      <tr>
        <td id="L599" class="blob-num js-line-number" data-line-number="599"></td>
        <td id="LC599" class="blob-code blob-code-inner js-file-line">				}</td>
      </tr>
      <tr>
        <td id="L600" class="blob-num js-line-number" data-line-number="600"></td>
        <td id="LC600" class="blob-code blob-code-inner js-file-line">			}</td>
      </tr>
      <tr>
        <td id="L601" class="blob-num js-line-number" data-line-number="601"></td>
        <td id="LC601" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L602" class="blob-num js-line-number" data-line-number="602"></td>
        <td id="LC602" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">return</span>(<span class="pl-k">list</span>(<span class="pl-smi">newColors</span>,<span class="pl-smi">colori</span>,unique(<span class="pl-smi">newColors</span>) ))</td>
      </tr>
      <tr>
        <td id="L603" class="blob-num js-line-number" data-line-number="603"></td>
        <td id="LC603" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L604" class="blob-num js-line-number" data-line-number="604"></td>
        <td id="LC604" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L605" class="blob-num js-line-number" data-line-number="605"></td>
        <td id="LC605" class="blob-code blob-code-inner js-file-line"><span class="pl-v">clusterTrackingPlot</span> <span class="pl-k">=</span> <span class="pl-k">function</span>(<span class="pl-smi">m</span>){</td>
      </tr>
      <tr>
        <td id="L606" class="blob-num js-line-number" data-line-number="606"></td>
        <td id="LC606" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>description: plots cluster tracking plot</span></td>
      </tr>
      <tr>
        <td id="L607" class="blob-num js-line-number" data-line-number="607"></td>
        <td id="LC607" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>input: m - matrix where rows are k, columns are samples, and values are cluster assignments.</span></td>
      </tr>
      <tr>
        <td id="L608" class="blob-num js-line-number" data-line-number="608"></td>
        <td id="LC608" class="blob-code blob-code-inner js-file-line">  plot(<span class="pl-c1">NULL</span>,<span class="pl-v">xlim</span><span class="pl-k">=</span>c(<span class="pl-k">-</span><span class="pl-c1">0.1</span>,<span class="pl-c1">1</span>),<span class="pl-v">ylim</span><span class="pl-k">=</span>c(<span class="pl-c1">0</span>,<span class="pl-c1">1</span>),<span class="pl-v">axes</span><span class="pl-k">=</span><span class="pl-c1">FALSE</span>,<span class="pl-v">xlab</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>samples<span class="pl-pds">&quot;</span></span>,<span class="pl-v">ylab</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>k<span class="pl-pds">&quot;</span></span>,<span class="pl-v">main</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>tracking plot<span class="pl-pds">&quot;</span></span>)</td>
      </tr>
      <tr>
        <td id="L609" class="blob-num js-line-number" data-line-number="609"></td>
        <td id="LC609" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">for</span>(<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span>nrow(<span class="pl-smi">m</span>)){</td>
      </tr>
      <tr>
        <td id="L610" class="blob-num js-line-number" data-line-number="610"></td>
        <td id="LC610" class="blob-code blob-code-inner js-file-line">    rect(  <span class="pl-v">xleft</span><span class="pl-k">=</span>seq(<span class="pl-c1">0</span>,<span class="pl-c1">1</span><span class="pl-k">-</span><span class="pl-c1">1</span><span class="pl-k">/</span>ncol(<span class="pl-smi">m</span>),<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span><span class="pl-k">/</span>ncol(<span class="pl-smi">m</span>)),  <span class="pl-v">ybottom</span><span class="pl-k">=</span>rep(<span class="pl-c1">1</span><span class="pl-k">-</span><span class="pl-smi">i</span><span class="pl-k">/</span>nrow(<span class="pl-smi">m</span>),ncol(<span class="pl-smi">m</span>)) , <span class="pl-v">xright</span><span class="pl-k">=</span>seq(<span class="pl-c1">1</span><span class="pl-k">/</span>ncol(<span class="pl-smi">m</span>),<span class="pl-c1">1</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span><span class="pl-k">/</span>ncol(<span class="pl-smi">m</span>)), <span class="pl-v">ytop</span><span class="pl-k">=</span>rep(<span class="pl-c1">1</span><span class="pl-k">-</span>(<span class="pl-smi">i</span><span class="pl-k">-</span><span class="pl-c1">1</span>)<span class="pl-k">/</span>nrow(<span class="pl-smi">m</span>),ncol(<span class="pl-smi">m</span>)), <span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">m</span>[<span class="pl-smi">i</span>,],<span class="pl-v">border</span><span class="pl-k">=</span><span class="pl-c1">NA</span>)   </td>
      </tr>
      <tr>
        <td id="L611" class="blob-num js-line-number" data-line-number="611"></td>
        <td id="LC611" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L612" class="blob-num js-line-number" data-line-number="612"></td>
        <td id="LC612" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>hatch lines to indicate samples</span></td>
      </tr>
      <tr>
        <td id="L613" class="blob-num js-line-number" data-line-number="613"></td>
        <td id="LC613" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">xl</span> <span class="pl-k">=</span> seq(<span class="pl-c1">0</span>,<span class="pl-c1">1</span><span class="pl-k">-</span><span class="pl-c1">1</span><span class="pl-k">/</span>ncol(<span class="pl-smi">m</span>),<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span><span class="pl-k">/</span>ncol(<span class="pl-smi">m</span>))</td>
      </tr>
      <tr>
        <td id="L614" class="blob-num js-line-number" data-line-number="614"></td>
        <td id="LC614" class="blob-code blob-code-inner js-file-line">  segments(  <span class="pl-smi">xl</span>, rep(<span class="pl-k">-</span><span class="pl-c1">0.1</span>,ncol(<span class="pl-smi">m</span>)) , <span class="pl-smi">xl</span>, rep(<span class="pl-c1">0</span>,ncol(<span class="pl-smi">m</span>)), <span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-s"><span class="pl-pds">&quot;</span>black<span class="pl-pds">&quot;</span></span>)    <span class="pl-c"><span class="pl-c">#</span>** alt white and black color?</span></td>
      </tr>
      <tr>
        <td id="L615" class="blob-num js-line-number" data-line-number="615"></td>
        <td id="LC615" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">ypos</span> <span class="pl-k">=</span> seq(<span class="pl-c1">1</span>,<span class="pl-c1">0</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-k">-</span><span class="pl-c1">1</span><span class="pl-k">/</span>nrow(<span class="pl-smi">m</span>))<span class="pl-k">-</span><span class="pl-c1">1</span><span class="pl-k">/</span>(<span class="pl-c1">2</span><span class="pl-k">*</span>nrow(<span class="pl-smi">m</span>))</td>
      </tr>
      <tr>
        <td id="L616" class="blob-num js-line-number" data-line-number="616"></td>
        <td id="LC616" class="blob-code blob-code-inner js-file-line">  text(<span class="pl-v">x</span><span class="pl-k">=</span><span class="pl-k">-</span><span class="pl-c1">0.1</span>,<span class="pl-v">y</span><span class="pl-k">=</span><span class="pl-smi">ypos</span>[<span class="pl-k">-</span>length(<span class="pl-smi">ypos</span>)],<span class="pl-v">labels</span><span class="pl-k">=</span>seq(<span class="pl-c1">2</span>,nrow(<span class="pl-smi">m</span>)<span class="pl-k">+</span><span class="pl-c1">1</span>,<span class="pl-v">by</span><span class="pl-k">=</span><span class="pl-c1">1</span>))</td>
      </tr>
      <tr>
        <td id="L617" class="blob-num js-line-number" data-line-number="617"></td>
        <td id="LC617" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L618" class="blob-num js-line-number" data-line-number="618"></td>
        <td id="LC618" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L619" class="blob-num js-line-number" data-line-number="619"></td>
        <td id="LC619" class="blob-code blob-code-inner js-file-line"><span class="pl-v">triangle</span> <span class="pl-k">=</span> <span class="pl-k">function</span>(<span class="pl-smi">m</span>,<span class="pl-v">mode</span><span class="pl-k">=</span><span class="pl-c1">1</span>){</td>
      </tr>
      <tr>
        <td id="L620" class="blob-num js-line-number" data-line-number="620"></td>
        <td id="LC620" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>mode=1 for CDF, vector of lower triangle.</span></td>
      </tr>
      <tr>
        <td id="L621" class="blob-num js-line-number" data-line-number="621"></td>
        <td id="LC621" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>mode==3 for full matrix.</span></td>
      </tr>
      <tr>
        <td id="L622" class="blob-num js-line-number" data-line-number="622"></td>
        <td id="LC622" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>mode==2 for calcICL; nonredundant half matrix coun</span></td>
      </tr>
      <tr>
        <td id="L623" class="blob-num js-line-number" data-line-number="623"></td>
        <td id="LC623" class="blob-code blob-code-inner js-file-line">  <span class="pl-c"><span class="pl-c">#</span>mode!=1 for summary </span></td>
      </tr>
      <tr>
        <td id="L624" class="blob-num js-line-number" data-line-number="624"></td>
        <td id="LC624" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">n</span><span class="pl-k">=</span>dim(<span class="pl-smi">m</span>)[<span class="pl-c1">1</span>]</td>
      </tr>
      <tr>
        <td id="L625" class="blob-num js-line-number" data-line-number="625"></td>
        <td id="LC625" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">nm</span> <span class="pl-k">=</span> <span class="pl-k">matrix</span>(<span class="pl-c1">0</span>,<span class="pl-v">ncol</span><span class="pl-k">=</span><span class="pl-smi">n</span>,<span class="pl-v">nrow</span><span class="pl-k">=</span><span class="pl-smi">n</span>)</td>
      </tr>
      <tr>
        <td id="L626" class="blob-num js-line-number" data-line-number="626"></td>
        <td id="LC626" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">fm</span> <span class="pl-k">=</span> <span class="pl-smi">m</span></td>
      </tr>
      <tr>
        <td id="L627" class="blob-num js-line-number" data-line-number="627"></td>
        <td id="LC627" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L628" class="blob-num js-line-number" data-line-number="628"></td>
        <td id="LC628" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L629" class="blob-num js-line-number" data-line-number="629"></td>
        <td id="LC629" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">nm</span>[upper.tri(<span class="pl-smi">nm</span>)] <span class="pl-k">=</span> <span class="pl-smi">m</span>[upper.tri(<span class="pl-smi">m</span>)] <span class="pl-c"><span class="pl-c">#</span>only upper half</span></td>
      </tr>
      <tr>
        <td id="L630" class="blob-num js-line-number" data-line-number="630"></td>
        <td id="LC630" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L631" class="blob-num js-line-number" data-line-number="631"></td>
        <td id="LC631" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">fm</span> <span class="pl-k">=</span> t(<span class="pl-smi">nm</span>)<span class="pl-k">+</span><span class="pl-smi">nm</span></td>
      </tr>
      <tr>
        <td id="L632" class="blob-num js-line-number" data-line-number="632"></td>
        <td id="LC632" class="blob-code blob-code-inner js-file-line">  diag(<span class="pl-smi">fm</span>) <span class="pl-k">=</span> diag(<span class="pl-smi">m</span>)</td>
      </tr>
      <tr>
        <td id="L633" class="blob-num js-line-number" data-line-number="633"></td>
        <td id="LC633" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L634" class="blob-num js-line-number" data-line-number="634"></td>
        <td id="LC634" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">nm</span><span class="pl-k">=</span><span class="pl-smi">fm</span></td>
      </tr>
      <tr>
        <td id="L635" class="blob-num js-line-number" data-line-number="635"></td>
        <td id="LC635" class="blob-code blob-code-inner js-file-line">  <span class="pl-smi">nm</span>[upper.tri(<span class="pl-smi">nm</span>)] <span class="pl-k">=</span> <span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L636" class="blob-num js-line-number" data-line-number="636"></td>
        <td id="LC636" class="blob-code blob-code-inner js-file-line">  diag(<span class="pl-smi">nm</span>) <span class="pl-k">=</span> <span class="pl-c1">NA</span></td>
      </tr>
      <tr>
        <td id="L637" class="blob-num js-line-number" data-line-number="637"></td>
        <td id="LC637" class="blob-code blob-code-inner js-file-line">  <span class="pl-v">vm</span> <span class="pl-k">=</span> <span class="pl-smi">m</span>[lower.tri(<span class="pl-smi">nm</span>)]</td>
      </tr>
      <tr>
        <td id="L638" class="blob-num js-line-number" data-line-number="638"></td>
        <td id="LC638" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L639" class="blob-num js-line-number" data-line-number="639"></td>
        <td id="LC639" class="blob-code blob-code-inner js-file-line">  <span class="pl-k">if</span>(<span class="pl-smi">mode</span><span class="pl-k">==</span><span class="pl-c1">1</span>){</td>
      </tr>
      <tr>
        <td id="L640" class="blob-num js-line-number" data-line-number="640"></td>
        <td id="LC640" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">return</span>(<span class="pl-smi">vm</span>) <span class="pl-c"><span class="pl-c">#</span>vector 		</span></td>
      </tr>
      <tr>
        <td id="L641" class="blob-num js-line-number" data-line-number="641"></td>
        <td id="LC641" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">mode</span><span class="pl-k">==</span><span class="pl-c1">3</span>){</td>
      </tr>
      <tr>
        <td id="L642" class="blob-num js-line-number" data-line-number="642"></td>
        <td id="LC642" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">return</span>(<span class="pl-smi">fm</span>) <span class="pl-c"><span class="pl-c">#</span>return full matrix</span></td>
      </tr>
      <tr>
        <td id="L643" class="blob-num js-line-number" data-line-number="643"></td>
        <td id="LC643" class="blob-code blob-code-inner js-file-line">  }<span class="pl-k">else</span> <span class="pl-k">if</span>(<span class="pl-smi">mode</span> <span class="pl-k">==</span> <span class="pl-c1">2</span>){</td>
      </tr>
      <tr>
        <td id="L644" class="blob-num js-line-number" data-line-number="644"></td>
        <td id="LC644" class="blob-code blob-code-inner js-file-line">    <span class="pl-k">return</span>(<span class="pl-smi">nm</span>) <span class="pl-c"><span class="pl-c">#</span>returns lower triangle and no diagonal. no double counts.</span></td>
      </tr>
      <tr>
        <td id="L645" class="blob-num js-line-number" data-line-number="645"></td>
        <td id="LC645" class="blob-code blob-code-inner js-file-line">  }</td>
      </tr>
      <tr>
        <td id="L646" class="blob-num js-line-number" data-line-number="646"></td>
        <td id="LC646" class="blob-code blob-code-inner js-file-line">  </td>
      </tr>
      <tr>
        <td id="L647" class="blob-num js-line-number" data-line-number="647"></td>
        <td id="LC647" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
      <tr>
        <td id="L648" class="blob-num js-line-number" data-line-number="648"></td>
        <td id="LC648" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L649" class="blob-num js-line-number" data-line-number="649"></td>
        <td id="LC649" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L650" class="blob-num js-line-number" data-line-number="650"></td>
        <td id="LC650" class="blob-code blob-code-inner js-file-line"><span class="pl-v">rankedBarPlot</span><span class="pl-k">=</span><span class="pl-k">function</span>(<span class="pl-smi">d</span>,<span class="pl-smi">myc</span>,<span class="pl-smi">cc</span>,<span class="pl-smi">title</span>){</td>
      </tr>
      <tr>
        <td id="L651" class="blob-num js-line-number" data-line-number="651"></td>
        <td id="LC651" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">colors</span> <span class="pl-k">=</span> rbind() <span class="pl-c"><span class="pl-c">#</span>each row is a barplot series</span></td>
      </tr>
      <tr>
        <td id="L652" class="blob-num js-line-number" data-line-number="652"></td>
        <td id="LC652" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">byRank</span> <span class="pl-k">=</span> cbind()</td>
      </tr>
      <tr>
        <td id="L653" class="blob-num js-line-number" data-line-number="653"></td>
        <td id="LC653" class="blob-code blob-code-inner js-file-line">
</td>
      </tr>
      <tr>
        <td id="L654" class="blob-num js-line-number" data-line-number="654"></td>
        <td id="LC654" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">spaceh</span> <span class="pl-k">=</span> <span class="pl-c1">0.1</span> <span class="pl-c"><span class="pl-c">#</span>space between bars</span></td>
      </tr>
      <tr>
        <td id="L655" class="blob-num js-line-number" data-line-number="655"></td>
        <td id="LC655" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">for</span>(<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">1</span><span class="pl-k">:</span>ncol(<span class="pl-smi">d</span>)){</td>
      </tr>
      <tr>
        <td id="L656" class="blob-num js-line-number" data-line-number="656"></td>
        <td id="LC656" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">byRank</span> <span class="pl-k">=</span> cbind(<span class="pl-smi">byRank</span>,sort(<span class="pl-smi">d</span>[,<span class="pl-smi">i</span>],<span class="pl-v">na.last</span><span class="pl-k">=</span><span class="pl-c1">F</span>))</td>
      </tr>
      <tr>
        <td id="L657" class="blob-num js-line-number" data-line-number="657"></td>
        <td id="LC657" class="blob-code blob-code-inner js-file-line">	  <span class="pl-v">colors</span> <span class="pl-k">=</span> rbind(<span class="pl-smi">colors</span>,order(<span class="pl-smi">d</span>[,<span class="pl-smi">i</span>],<span class="pl-v">na.last</span><span class="pl-k">=</span><span class="pl-c1">F</span>))</td>
      </tr>
      <tr>
        <td id="L658" class="blob-num js-line-number" data-line-number="658"></td>
        <td id="LC658" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L659" class="blob-num js-line-number" data-line-number="659"></td>
        <td id="LC659" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">maxH</span> <span class="pl-k">=</span> max(c(<span class="pl-c1">1.5</span>,apply(<span class="pl-smi">byRank</span>,<span class="pl-c1">2</span>,<span class="pl-smi">sum</span>)),<span class="pl-v">na.rm</span><span class="pl-k">=</span><span class="pl-c1">T</span>) <span class="pl-c"><span class="pl-c">#</span>maximum height of graph</span></td>
      </tr>
      <tr>
        <td id="L660" class="blob-num js-line-number" data-line-number="660"></td>
        <td id="LC660" class="blob-code blob-code-inner js-file-line">	</td>
      </tr>
      <tr>
        <td id="L661" class="blob-num js-line-number" data-line-number="661"></td>
        <td id="LC661" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">#</span>barplot largest to smallest so that smallest is in front.</span></td>
      </tr>
      <tr>
        <td id="L662" class="blob-num js-line-number" data-line-number="662"></td>
        <td id="LC662" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">barp</span> <span class="pl-k">=</span> barplot( apply(<span class="pl-smi">byRank</span>,<span class="pl-c1">2</span>,<span class="pl-smi">sum</span>) ,  <span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">myc</span>[<span class="pl-smi">colors</span>[,<span class="pl-c1">1</span>]] ,<span class="pl-v">space</span><span class="pl-k">=</span><span class="pl-smi">spaceh</span>,<span class="pl-v">ylim</span><span class="pl-k">=</span>c(<span class="pl-c1">0</span>,<span class="pl-smi">maxH</span>),<span class="pl-v">main</span><span class="pl-k">=</span>paste(<span class="pl-s"><span class="pl-pds">&quot;</span>item-consensus<span class="pl-pds">&quot;</span></span>, <span class="pl-smi">title</span>),<span class="pl-v">border</span><span class="pl-k">=</span><span class="pl-c1">NA</span>,<span class="pl-v">las</span><span class="pl-k">=</span><span class="pl-c1">1</span>  )</td>
      </tr>
      <tr>
        <td id="L663" class="blob-num js-line-number" data-line-number="663"></td>
        <td id="LC663" class="blob-code blob-code-inner js-file-line">	<span class="pl-k">for</span>(<span class="pl-smi">i</span> <span class="pl-k">in</span> <span class="pl-c1">2</span><span class="pl-k">:</span>nrow(<span class="pl-smi">byRank</span>)){</td>
      </tr>
      <tr>
        <td id="L664" class="blob-num js-line-number" data-line-number="664"></td>
        <td id="LC664" class="blob-code blob-code-inner js-file-line">	  barplot( apply(<span class="pl-k">matrix</span>(<span class="pl-smi">byRank</span>[<span class="pl-smi">i</span><span class="pl-k">:</span>nrow(<span class="pl-smi">byRank</span>),],<span class="pl-v">ncol</span><span class="pl-k">=</span>ncol(<span class="pl-smi">byRank</span>))  ,<span class="pl-c1">2</span>,<span class="pl-smi">sum</span>), <span class="pl-v">space</span><span class="pl-k">=</span><span class="pl-smi">spaceh</span>,<span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">myc</span>[<span class="pl-smi">colors</span>[,<span class="pl-smi">i</span>]],<span class="pl-v">ylim</span><span class="pl-k">=</span>c(<span class="pl-c1">0</span>,<span class="pl-smi">maxH</span>), <span class="pl-v">add</span><span class="pl-k">=</span><span class="pl-c1">T</span>,<span class="pl-v">border</span><span class="pl-k">=</span><span class="pl-c1">NA</span>,<span class="pl-v">las</span><span class="pl-k">=</span><span class="pl-c1">1</span>  )</td>
      </tr>
      <tr>
        <td id="L665" class="blob-num js-line-number" data-line-number="665"></td>
        <td id="LC665" class="blob-code blob-code-inner js-file-line">	}</td>
      </tr>
      <tr>
        <td id="L666" class="blob-num js-line-number" data-line-number="666"></td>
        <td id="LC666" class="blob-code blob-code-inner js-file-line">	<span class="pl-v">xr</span><span class="pl-k">=</span>seq(<span class="pl-smi">spaceh</span>,ncol(<span class="pl-smi">d</span>)<span class="pl-k">+</span>ncol(<span class="pl-smi">d</span>)<span class="pl-k">*</span><span class="pl-smi">spaceh</span>,(ncol(<span class="pl-smi">d</span>)<span class="pl-k">+</span>ncol(<span class="pl-smi">d</span>)<span class="pl-k">*</span><span class="pl-smi">spaceh</span>)<span class="pl-k">/</span>ncol(<span class="pl-smi">d</span>)  )</td>
      </tr>
      <tr>
        <td id="L667" class="blob-num js-line-number" data-line-number="667"></td>
        <td id="LC667" class="blob-code blob-code-inner js-file-line">	<span class="pl-c"><span class="pl-c">#</span>class labels as asterisks</span></td>
      </tr>
      <tr>
        <td id="L668" class="blob-num js-line-number" data-line-number="668"></td>
        <td id="LC668" class="blob-code blob-code-inner js-file-line">	text(<span class="pl-s"><span class="pl-pds">&quot;</span>*<span class="pl-pds">&quot;</span></span>,<span class="pl-v">x</span><span class="pl-k">=</span><span class="pl-smi">xr</span><span class="pl-k">+</span><span class="pl-c1">0.5</span>,<span class="pl-v">y</span><span class="pl-k">=</span><span class="pl-smi">maxH</span>,<span class="pl-v">col</span><span class="pl-k">=</span><span class="pl-smi">myc</span>[<span class="pl-smi">cc</span>],<span class="pl-v">cex</span><span class="pl-k">=</span><span class="pl-c1">1.4</span>) <span class="pl-c"><span class="pl-c">#</span>rect(xr,1.4,xr+1,1.5,col=myc[cc] )</span></td>
      </tr>
      <tr>
        <td id="L669" class="blob-num js-line-number" data-line-number="669"></td>
        <td id="LC669" class="blob-code blob-code-inner js-file-line">}</td>
      </tr>
</table>

  <details class="details-reset details-overlay BlobToolbar position-absolute js-file-line-actions dropdown d-none" aria-hidden="true">
    <summary class="btn-octicon ml-0 px-2 p-0 bg-white border border-gray-dark rounded-1" aria-label="Inline file action toolbar">
      <svg class="octicon octicon-kebab-horizontal" viewBox="0 0 13 16" version="1.1" width="13" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M1.5 9a1.5 1.5 0 1 0 0-3 1.5 1.5 0 0 0 0 3zm5 0a1.5 1.5 0 1 0 0-3 1.5 1.5 0 0 0 0 3zM13 7.5a1.5 1.5 0 1 1-3 0 1.5 1.5 0 0 1 3 0z"/></svg>
    </summary>
    <details-menu>
      <ul class="BlobToolbar-dropdown dropdown-menu dropdown-menu-se mt-2" style="width:185px">
        <li><clipboard-copy role="menuitem" class="dropdown-item" id="js-copy-lines" style="cursor:pointer;" data-original-text="Copy lines">Copy lines</clipboard-copy></li>
        <li><clipboard-copy role="menuitem" class="dropdown-item" id="js-copy-permalink" style="cursor:pointer;" data-original-text="Copy permalink">Copy permalink</clipboard-copy></li>
        <li><a class="dropdown-item js-update-url-with-hash" id="js-view-git-blame" role="menuitem" href="/renzhonglu/ConsensusClusterPlus/blame/08bb6aa1e93619807738d015d8c40f4bf301266b/R/ConsensusClusterPlus.R">View git blame</a></li>
      </ul>
    </details-menu>
  </details>

  </div>

    </div>

  

  <details class="details-reset details-overlay details-overlay-dark">
    <summary data-hotkey="l" aria-label="Jump to line"></summary>
    <details-dialog class="Box Box--overlay d-flex flex-column anim-fade-in fast linejump" aria-label="Jump to line">
      <!-- '"` --><!-- </textarea></xmp> --></option></form><form class="js-jump-to-line-form Box-body d-flex" action="" accept-charset="UTF-8" method="get"><input name="utf8" type="hidden" value="&#x2713;" />
        <input class="form-control flex-auto mr-3 linejump-input js-jump-to-line-field" type="text" placeholder="Jump to line&hellip;" aria-label="Jump to line" autofocus>
        <button type="submit" class="btn" data-close-dialog>Go</button>
</form>    </details-dialog>
  </details>



  </div>
  <div class="modal-backdrop js-touch-events"></div>
</div>

    </main>
  </div>
  

  </div>

        
<div class="footer container-lg width-full px-3" role="contentinfo">
  <div class="position-relative d-flex flex-justify-between pt-6 pb-2 mt-6 f6 text-gray border-top border-gray-light ">
    <ul class="list-style-none d-flex flex-wrap ">
      <li class="mr-3">&copy; 2019 <span title="0.16575s from unicorn-55c9fbbdc5-dbhvq">GitHub</span>, Inc.</li>
        <li class="mr-3"><a data-ga-click="Footer, go to terms, text:terms" href="https://github.com/site/terms">Terms</a></li>
        <li class="mr-3"><a data-ga-click="Footer, go to privacy, text:privacy" href="https://github.com/site/privacy">Privacy</a></li>
        <li class="mr-3"><a data-ga-click="Footer, go to security, text:security" href="https://github.com/security">Security</a></li>
        <li class="mr-3"><a href="https://githubstatus.com/" data-ga-click="Footer, go to status, text:status">Status</a></li>
        <li><a data-ga-click="Footer, go to help, text:help" href="https://help.github.com">Help</a></li>
    </ul>

    <a aria-label="Homepage" title="GitHub" class="footer-octicon mx-lg-4" href="https://github.com">
      <svg height="24" class="octicon octicon-mark-github" viewBox="0 0 16 16" version="1.1" width="24" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"/></svg>
</a>
   <ul class="list-style-none d-flex flex-wrap ">
        <li class="mr-3"><a data-ga-click="Footer, go to contact, text:contact" href="https://github.com/contact">Contact GitHub</a></li>
        <li class="mr-3"><a href="https://github.com/pricing" data-ga-click="Footer, go to Pricing, text:Pricing">Pricing</a></li>
      <li class="mr-3"><a href="https://developer.github.com" data-ga-click="Footer, go to api, text:api">API</a></li>
      <li class="mr-3"><a href="https://training.github.com" data-ga-click="Footer, go to training, text:training">Training</a></li>
        <li class="mr-3"><a href="https://github.blog" data-ga-click="Footer, go to blog, text:blog">Blog</a></li>
        <li><a data-ga-click="Footer, go to about, text:about" href="https://github.com/about">About</a></li>

    </ul>
  </div>
  <div class="d-flex flex-justify-center pb-6">
    <span class="f6 text-gray-light"></span>
  </div>
</div>



  <div id="ajax-error-message" class="ajax-error-message flash flash-error">
    <svg class="octicon octicon-alert" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"/></svg>
    <button type="button" class="flash-close js-ajax-error-dismiss" aria-label="Dismiss error">
      <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
    </button>
    You can’t perform that action at this time.
  </div>


    <script crossorigin="anonymous" integrity="sha512-5XCib/MAIgqVSU6/fcOJEezL7AmEQhcmcvMnO/uRQLBGEPQm2jFN2QHDw+49RhVBHdjYAQtBq+aKcsusMsU2Kw==" type="application/javascript" src="https://github.githubassets.com/assets/compat-bootstrap-9148e4a7.js"></script>
    <script crossorigin="anonymous" integrity="sha512-IZCdsT3BCdoe3zLHNBvhGH1FjbGHtj+WInKlc6pqAwOglygvSKUfOPyxxtJco1AAVsVk1ldvDmCJyuIYNZhy7g==" type="application/javascript" src="https://github.githubassets.com/assets/frameworks-b805de66.js"></script>
    
    <script crossorigin="anonymous" async="async" integrity="sha512-ScnmScO/jzVnXySogyQfC0UfuKIFipm/GKpRDt4RKKbBxPR3jAkAULvg73h21bbY7MUZyhP5Ti6Amcbd9fH4Gw==" type="application/javascript" src="https://github.githubassets.com/assets/github-bootstrap-0179ec3c.js"></script>
    
    
    
  <div class="js-stale-session-flash stale-session-flash flash flash-warn flash-banner d-none">
    <svg class="octicon octicon-alert" viewBox="0 0 16 16" version="1.1" width="16" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M8.893 1.5c-.183-.31-.52-.5-.887-.5s-.703.19-.886.5L.138 13.499a.98.98 0 0 0 0 1.001c.193.31.53.501.886.501h13.964c.367 0 .704-.19.877-.5a1.03 1.03 0 0 0 .01-1.002L8.893 1.5zm.133 11.497H6.987v-2.003h2.039v2.003zm0-3.004H6.987V5.987h2.039v4.006z"/></svg>
    <span class="signed-in-tab-flash">You signed in with another tab or window. <a href="">Reload</a> to refresh your session.</span>
    <span class="signed-out-tab-flash">You signed out in another tab or window. <a href="">Reload</a> to refresh your session.</span>
  </div>
  <template id="site-details-dialog">
  <details class="details-reset details-overlay details-overlay-dark lh-default text-gray-dark" open>
    <summary aria-haspopup="dialog" aria-label="Close dialog"></summary>
    <details-dialog class="Box Box--overlay d-flex flex-column anim-fade-in fast">
      <button class="Box-btn-octicon m-0 btn-octicon position-absolute right-0 top-0" type="button" aria-label="Close dialog" data-close-dialog>
        <svg class="octicon octicon-x" viewBox="0 0 12 16" version="1.1" width="12" height="16" aria-hidden="true"><path fill-rule="evenodd" d="M7.48 8l3.75 3.75-1.48 1.48L6 9.48l-3.75 3.75-1.48-1.48L4.52 8 .77 4.25l1.48-1.48L6 6.52l3.75-3.75 1.48 1.48L7.48 8z"/></svg>
      </button>
      <div class="octocat-spinner my-6 js-details-dialog-spinner"></div>
    </details-dialog>
  </details>
</template>

  <div class="Popover js-hovercard-content position-absolute" style="display: none; outline: none;" tabindex="0">
  <div class="Popover-message Popover-message--bottom-left Popover-message--large Box box-shadow-large" style="width:360px;">
  </div>
</div>

  <div aria-live="polite" class="js-global-screen-reader-notice sr-only"></div>

  </body>
</html>

