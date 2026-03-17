# ChemBind — Deployment Guide

## Prerequisites

- Firebase project with Firestore + Authentication (Google provider) enabled
- Firebase service account JSON (download from Firebase Console → Project Settings → Service Accounts)
- Git repository pushed to GitHub

---

## 1. Managed Redis (Upstash)

1. Go to [Upstash Console](https://console.upstash.com/)
2. Create a new Redis database (free tier)
3. Copy the **`rediss://...`** connection URL (TLS-enabled)
4. You'll use this as `REDIS_URL` for both API and Worker services

---

## 2. Deploy API Service

### Option A: Railway

1. [railway.app](https://railway.app/) → New Project → Deploy from GitHub repo
2. Set root directory: `/` (Dockerfile is at `backend/Dockerfile` and uses context `.`)
3. Railway auto-detects `backend/railway.json`
4. Add environment variables:

| Variable | Value |
|---|---|
| `FIREBASE_SERVICE_ACCOUNT_JSON` | Paste entire service-account JSON (single line) |
| `REDIS_URL` | `rediss://default:PASSWORD@HOST:PORT` from Upstash |
| `CELERY_BROKER_URL` | Same as `REDIS_URL` |
| `CELERY_RESULT_BACKEND` | Same as `REDIS_URL` (Celery auto-appends `/1`) |
| `RATE_LIMIT_STORAGE_URI` | Same as `REDIS_URL` (or leave empty for in-memory) |
| `CORS_ORIGINS` | `https://your-app.vercel.app` (update after frontend deploy) |
| `TRUSTED_HOSTS` | `your-api.up.railway.app,localhost` |
| `SENTRY_DSN` | Your Sentry DSN (or leave empty) |
| `ENVIRONMENT` | `production` |
| `DEBUG` | `0` |

5. Deploy. Verify: `curl https://YOUR-API.up.railway.app/api/health` → `{"ok": true}`

### Option B: Render

1. [render.com](https://render.com/) → New Web Service → Connect GitHub repo
2. Root directory: leave empty (uses repo root)
3. Environment: Docker
4. Dockerfile path: `backend/Dockerfile`
5. Docker context: `.`
6. Health check path: `/api/health`
7. Add the same environment variables as Railway (above)
8. Deploy. Verify health endpoint.

---

## 3. Deploy Worker

### Railway
1. Same project → Add new service → Same GitHub repo
2. Override start command: `bash -lc "chmod +x /app/worker.sh && /app/worker.sh"`
3. Same env vars as API (especially `REDIS_URL`, `FIREBASE_SERVICE_ACCOUNT_JSON`)
4. No port needed (background worker)

### Render
The `render.yaml` blueprint already defines the worker service. If using the blueprint:
1. Connect repo → Render auto-discovers `render.yaml`
2. Fill in env var values in the Render dashboard

If creating manually:
1. New Background Worker → Docker → same Dockerfile
2. Start command: `bash -lc "chmod +x /app/worker.sh && /app/worker.sh"`
3. Same env vars as API

---

## 4. Verify

```bash
# API health
curl -s https://YOUR-API/api/health
# Expected: {"ok":true}

# Check worker logs in Railway/Render dashboard
# Expected: "ready" message, no Redis connection errors
```

---

## Environment Variables Reference

See `backend/.env.example` for all backend variables with descriptions.
See `frontend/.env.example` for all frontend variables.
