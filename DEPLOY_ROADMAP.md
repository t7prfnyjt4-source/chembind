# Segments 49-51: Deploy Roadmap Features + Regression + QA

## Segment 49 — Deploy Roadmap Features

### Backend Deployment
1. Rebuild Docker image (includes Vina, Meeko, ProLIF):
   ```bash
   docker build -t chembind-backend ./backend
   ```

2. Set feature flags in production environment:
   ```
   ENABLE_SIMILARITY_SEARCH=true
   ENABLE_CONFORMERS=true
   ENABLE_DOCKING=true
   ENABLE_EXPORT=true
   ENABLE_ANNOTATIONS=true
   ```

3. Deploy backend service (Railway/Render)
4. Deploy worker service with same image
5. Deploy updated Firestore rules:
   ```bash
   firebase deploy --only firestore:rules
   ```

### Frontend Deployment
1. Deploy to Vercel:
   ```bash
   cd frontend && npx vercel --prod
   ```
2. Verify VITE_API_BASE points to production API

### Post-Deploy Checks
- [ ] `curl https://YOUR-API/api/health` → 200
- [ ] Worker logs show "ready"
- [ ] No new Sentry errors

---

## Segment 50 — Full Regression Test Checklist

### CORE
- [ ] GET /api/health → 200
- [ ] POST /api/analyze with valid SMILES → descriptors
- [ ] POST /api/batch with 5 rows → job completes
- [ ] Idempotency: same CSV re-upload returns same jobId
- [ ] Auth: 401 without token on protected endpoints
- [ ] Rate limit: 65× /api/analyze → 429
- [ ] Headers: X-Request-Id, HSTS, X-Frame-Options present

### SEARCH
- [ ] POST /api/similarity/search → ranked results
- [ ] POST /api/search mode=substructure → matches
- [ ] Empty history → empty results

### CONFORMERS
- [ ] POST /api/conformers → conformers with coords
- [ ] Second call → cached=true
- [ ] >80 atoms → error

### DOCKING
- [ ] POST /api/docking/jobs with PDB + SMILES → job created
- [ ] Job progresses: queued → preparing → docking → completed
- [ ] GET /api/docking/jobs/{id}/poses → sorted poses
- [ ] Non-enterprise user → 402
- [ ] Interactions stored per pose

### EXPORT
- [ ] GET /api/export?format=mol → downloads .mol
- [ ] GET /api/export?format=sdf → downloads .sdf
- [ ] GET /api/export?format=cdxml → downloads .cdxml
- [ ] Non-professional user → 402

### GALLERY
- [ ] Batch items have thumbnail field
- [ ] Gallery grid renders
- [ ] Click-through works

### COLLAB
- [ ] POST /api/annotations → creates annotation
- [ ] GET /api/annotations → lists user's annotations
- [ ] PUT /api/annotations/{id} → updates
- [ ] POST /api/annotations/{id}/share → returns shareId
- [ ] GET /api/shared/{shareId} → public view (no auth)

### TIERS
- [ ] basic user blocked from enterprise features → 402
- [ ] 402 messages include current and required tier

---

## Segment 51 — Docking + Export QA

### Docking QA (5 known protein-ligand pairs)
Test with these well-known systems:
1. 1iep + Imatinib (Gleevec)
2. 4ey7 + Erlotinib
3. 1err + Raloxifene
4. 3eml + Sorafenib
5. 2rh1 + Carazolol

Verify: poses have negative scores, binding makes chemical sense.

### Export QA
Generate CDXML for 20 molecules, open in ChemDraw to verify rendering.

### Gallery QA
Run batch of 50 diverse molecules, verify thumbnails and scroll.
