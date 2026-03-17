import React, { useEffect, useMemo, useState } from "react";
import {
  collection,
  limit,
  onSnapshot,
  orderBy,
  query,
  Timestamp,
} from "firebase/firestore";

import { db, auth } from "../firebase/config";

type AnalysisDoc = {
  smiles?: string;
  canonical?: string;
  valid?: boolean;
  mw?: number | null;
  logp?: number | null;
  createdAt?: Timestamp | { seconds: number } | null;
};

function fmtNum(x: unknown, digits = 2): string {
  if (x === null || x === undefined) return "—";
  const n = typeof x === "number" ? x : Number(x);
  if (!Number.isFinite(n)) return "—";
  return n.toFixed(digits);
}

function tsToDate(x: AnalysisDoc["createdAt"]): Date | null {
  if (!x) return null;
  if (x instanceof Timestamp) return x.toDate();
  if (typeof x === "object" && x !== null && "seconds" in x) {
    const sec = (x as { seconds: number }).seconds;
    if (Number.isFinite(sec)) return new Date(sec * 1000);
  }
  return null;
}

function Spinner({ label = "Loading…" }: { label?: string }) {
  return (
    <div style={{ display: "flex", gap: 8, alignItems: "center", fontSize: 14 }}>
      <span
        style={{
          width: 14,
          height: 14,
          border: "2px solid rgba(0,0,0,0.35)",
          borderTopColor: "transparent",
          borderRadius: "50%",
          display: "inline-block",
          animation: "spin 0.9s linear infinite",
        }}
      />
      <span>{label}</span>
      <style>
        {`
          @keyframes spin {
            from { transform: rotate(0deg); }
            to { transform: rotate(360deg); }
          }
        `}
      </style>
    </div>
  );
}

function Badge({ ok }: { ok: boolean }) {
  return (
    <span
      style={{
        fontSize: 12,
        padding: "2px 8px",
        borderRadius: 999,
        background: ok ? "#d1fae5" : "#fee2e2",
        color: ok ? "#065f46" : "#991b1b",
        fontWeight: 600,
      }}
    >
      {ok ? "valid" : "invalid"}
    </span>
  );
}

export default function AnalysisHistory() {
  const user = auth.currentUser;

  // Start true; turn false inside snapshot callback -> avoids setState-in-effect
  const [loading, setLoading] = useState(true);
  const [loadedOnce, setLoadedOnce] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [docs, setDocs] = useState<Array<{ id: string; data: AnalysisDoc }>>(
    []
  );
  const [expandedId, setExpandedId] = useState<string | null>(null);

  const emptyState = useMemo(() => {
    if (!user) {
      return { title: "Sign in to save", hint: "Log in to see your analysis history." };
    }
    if (loadedOnce && !loading && docs.length === 0) {
      return {
        title: "Analyze your first molecule",
        hint: "Run an analysis and it will appear here in real time.",
      };
    }
    return null;
  }, [user, loadedOnce, loading, docs.length]);

  useEffect(() => {
    if (!user) return;

    const colRef = collection(db, `users/${user.uid}/analyses`);
    const q = query(colRef, orderBy("createdAt", "desc"), limit(50));

    const unsub = onSnapshot(
      q,
      (snap) => {
        const results = snap.docs.map((d) => ({
          id: d.id,
          data: d.data() as AnalysisDoc,
        }));
        setDocs(results);
        setError(null);
        setLoading(false);
        setLoadedOnce(true);
      },
      (e: unknown) => {
        const msg = e instanceof Error ? e.message : "Failed to load history";
        setError(msg);
        setLoading(false);
        setLoadedOnce(true);
      }
    );

    return () => unsub();
  }, [user]);

  if (!user) {
    return (
      <div style={{ padding: 12 }}>
        <h3 style={{ margin: "0 0 6px 0" }}>{emptyState?.title}</h3>
        <p style={{ margin: 0, opacity: 0.8 }}>{emptyState?.hint}</p>
      </div>
    );
  }

  return (
    <div style={{ padding: 12 }}>
      <div
        style={{
          display: "flex",
          alignItems: "center",
          justifyContent: "space-between",
          gap: 12,
          marginBottom: 12,
        }}
      >
        <h2 style={{ margin: 0 }}>History</h2>
        {loading ? <Spinner /> : null}
      </div>

      {error ? (
        <div
          style={{
            border: "1px solid #fecaca",
            background: "#fef2f2",
            color: "#991b1b",
            padding: 10,
            borderRadius: 8,
            marginBottom: 12,
            fontSize: 14,
          }}
        >
          {error}
        </div>
      ) : null}

      {emptyState ? (
        <div
          style={{
            border: "1px solid #ddd",
            padding: 12,
            borderRadius: 8,
            marginBottom: 12,
          }}
        >
          <div style={{ fontWeight: 700 }}>{emptyState.title}</div>
          <div style={{ opacity: 0.8, marginTop: 4 }}>{emptyState.hint}</div>
        </div>
      ) : null}

      <div style={{ display: "grid", gap: 12 }}>
        {docs.map(({ id, data }) => {
          const expanded = expandedId === id;
          const created = tsToDate(data.createdAt);
          const ok = !!data.valid;

          return (
            <div
              key={id}
              onClick={() => setExpandedId(expanded ? null : id)}
              style={{
                border: "1px solid #444",
                borderRadius: 8,
                padding: 12,
                cursor: "pointer",
              }}
            >
              <div
                style={{
                  display: "flex",
                  justifyContent: "space-between",
                  gap: 12,
                  alignItems: "flex-start",
                }}
              >
                <div style={{ minWidth: 0 }}>
                  <div style={{ fontSize: 13, opacity: 0.7 }}>SMILES</div>
                  <div style={{ fontFamily: "monospace", wordBreak: "break-all" }}>
                    {data.smiles ?? "—"}
                  </div>
                </div>
                <Badge ok={ok} />
              </div>

              <div style={{ marginTop: 8, fontSize: 14 }}>
                <strong>MW:</strong> {fmtNum(data.mw)}{" "}
                <span style={{ opacity: 0.6 }}>|</span>{" "}
                <strong>LogP:</strong> {fmtNum(data.logp)}
              </div>

              {expanded ? (
                <div style={{ marginTop: 10, borderTop: "1px solid #555", paddingTop: 10 }}>
                  <div style={{ fontSize: 14 }}>
                    <strong>Canonical:</strong>{" "}
                    <span style={{ fontFamily: "monospace", wordBreak: "break-all" }}>
                      {data.canonical ?? "—"}
                    </span>
                  </div>
                  <div style={{ fontSize: 14, marginTop: 6 }}>
                    <strong>Created:</strong> {created ? created.toLocaleString() : "—"}
                  </div>
                  <div style={{ fontSize: 12, opacity: 0.7, marginTop: 8 }}>
                    Click to collapse
                  </div>
                </div>
              ) : (
                <div style={{ fontSize: 12, opacity: 0.7, marginTop: 8 }}>
                  Click to expand
                </div>
              )}
            </div>
          );
        })}
      </div>
    </div>
  );
}
