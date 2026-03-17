# ChemBind — E2E Production Checklist

Run this checklist against your production deployment after Segments 19–20 are live.

## AUTH
- [ ] Login with Google → success, token received
- [ ] Logout → session cleared
- [ ] Refresh page → stays logged in (token refresh)
- [ ] Request without token → 401

## ANALYZE
- [ ] Valid SMILES (e.g., `CCO`) → descriptors returned
- [ ] Invalid SMILES (e.g., `INVALID`) → 400 error
- [ ] Salt SMILES (e.g., `[Na+].[Cl-]`) → descriptors returned
- [ ] Long SMILES (400 chars) → success
- [ ] Very long SMILES (>500 chars) → 400 error

## BATCH
- [ ] Upload CSV with 5 rows → job created, progress updates
- [ ] Job completes → all items have results
- [ ] Idempotency: re-submit same CSV → returns same job ID
- [ ] >500 rows → 413 or 400 error

## RATE LIMITING
- [ ] 65× rapid analyze requests → 429 after limit hit

## HEADERS
- [ ] Response has `X-Content-Type-Options: nosniff`
- [ ] Response has `X-Frame-Options: DENY`
- [ ] Response has `X-Request-Id`
- [ ] HTTPS response has `Strict-Transport-Security`

## WORKER
- [ ] Kill worker process → restart → in-progress job resumes/completes

## FIRESTORE RULES
- [ ] User A cannot read User B's analyses (test with different tokens)

## BACKUPS
- [ ] Firestore backup schedule is active (Firebase Console → Firestore → Backups)

## SENTRY
- [ ] Trigger a test error → verify it appears in Sentry dashboard

## UPTIME
- [ ] Set up UptimeRobot (or similar) on `https://YOUR-API/api/health` with 5-min interval
- [ ] Verify first ping succeeds

---

If all checks pass: **production is live**. Fix anything that fails before proceeding.
