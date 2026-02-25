# frontend/Dockerfile

FROM node:20-alpine AS build
WORKDIR /app

# copy only frontend manifests first (better cache)
COPY frontend/package*.json ./
RUN npm ci

# copy the actual frontend source
COPY frontend ./
RUN npm run build

FROM nginx:alpine
COPY --from=build /app/dist /usr/share/nginx/html

# SPA routing
RUN rm /etc/nginx/conf.d/default.conf
COPY frontend/nginx.conf /etc/nginx/conf.d/default.conf

EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
