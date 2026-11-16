# nicohadas.github.io

A minimal Jekyll site.

- `index.html`: home page (name, links, email)
- `about.md`: About page
- `_data/publications.yml`: publications, which `publications.html` renders
- `_posts/`: blog posts, named `YYYY-MM-DD-title.md`. Drafts go in `_drafts/`.
- `assets/css/style.css`: all styling

Local preview: `bundle install && bundle exec jekyll serve --drafts`

Pushing to `main` runs `.github/workflows/deploy.yml`, which builds the site into the `gh-pages` branch.
