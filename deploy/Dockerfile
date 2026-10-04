# Build context is the repository root: `fly deploy` from there.
FROM nginx:alpine
COPY deploy/nginx.conf /etc/nginx/nginx.conf
COPY web/index.html web/explorer.js web/census-map.js /usr/share/nginx/html/
COPY web/vendor/ /usr/share/nginx/html/vendor/
# Precompress once at build time; nginx serves the .gz with gzip_static.
RUN cd /usr/share/nginx/html && gzip -k9 explorer.js census-map.js index.html \
      vendor/katex/katex.min.js vendor/katex/katex.min.css vendor/katex/contrib/auto-render.min.js
EXPOSE 8080
CMD ["nginx", "-g", "daemon off;"]
