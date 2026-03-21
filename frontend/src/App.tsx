import { useEffect, useState } from "react";
import BatchUpload from "./components/BatchUpload";

import {
  type User,
  onAuthStateChanged,
  GoogleAuthProvider,
  signInWithPopup,
  signOut,
} from "firebase/auth";

import { auth } from "./firebase/config";
import AnalysisHistory from "./components/AnalysisHistory";
import MoleculeSearch from "./components/MoleculeSearch";
import ConformerPlayer from "./components/ConformerPlayer";

type Page = "analyze" | "batch" | "history" | "search";

function Spinner({ label = "Loading…" }: { label?: string }) {
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 8, fontSize: 14, opacity: 0.85 }}>
      <span
        style={{
          width: 14,
          height: 14,
          borderRadius: "50%",
          border: "2px solid rgba(0,0,0,0.25)",
          borderTopColor: "transparent",
          display: "inline-block",
          animation: "spin 0.8s linear infinite",
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

export default function App() {
  const [user, setUser] = useState<User | null>(null);
  const [authLoading, setAuthLoading] = useState(true);
  const [page, setPage] = useState<Page>("analyze");
  const [analyzeSmiles, setAnalyzeSmiles] = useState("");

  useEffect(() => {
    const unsub = onAuthStateChanged(auth, (u) => {
      setUser(u);
      setAuthLoading(false);
    });
    return () => unsub();
  }, []);

  const login = async () => {
    const provider = new GoogleAuthProvider();
    await signInWithPopup(auth, provider);
  };

  const logout = async () => {
    await signOut(auth);
  };

  if (authLoading) {
    return (
      <div style={{ padding: 24, fontFamily: "system-ui" }}>
        <Spinner label="Checking login…" />
      </div>
    );
  }

  return (
    <div style={{ padding: 24, fontFamily: "system-ui", maxWidth: 1100, margin: "0 auto" }}>
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          gap: 12,
          flexWrap: "wrap",
        }}
      >
        <h1 style={{ margin: 0 }}>Chembind</h1>

        <div style={{ display: "flex", gap: 8, alignItems: "center", flexWrap: "wrap" }}>
          {user ? (
            <>
              <span style={{ opacity: 0.8, fontSize: 14 }}>
                {user.email ?? user.uid}
              </span>
              <button onClick={logout}>Logout</button>
            </>
          ) : (
            <button onClick={login}>Login</button>
          )}
        </div>
      </div>

      <div style={{ display: "flex", gap: 8, marginTop: 16, flexWrap: "wrap" }}>
        <button onClick={() => setPage("analyze")}>Analyze</button>
        <button onClick={() => setPage("batch")}>Batch</button>
        <button onClick={() => setPage("history")}>History</button>
        <button onClick={() => setPage("search")}>Search</button>
      </div>

      <div
        style={{
          marginTop: 16,
          padding: 12,
          border: "1px solid #333",
          borderRadius: 8,
        }}
      >
        {page === "analyze" && (
          <div>
            <div style={{ fontWeight: 600, marginBottom: 6 }}>Analyze</div>
            <div style={{ display: "flex", gap: 8, marginBottom: 8 }}>
              <input
                type="text"
                value={analyzeSmiles}
                onChange={(e) => setAnalyzeSmiles(e.target.value)}
                placeholder="e.g. CCO"
                style={{ padding: "6px 10px", fontFamily: "monospace", flex: 1, maxWidth: 400 }}
              />
            </div>
            {analyzeSmiles.trim() && <ConformerPlayer smiles={analyzeSmiles.trim()} />}
          </div>
        )}

        {page === "batch" && (
          <div>
            <div style={{ fontWeight: 600, marginBottom: 6 }}>Batch</div>
            <BatchUpload token="" />
          </div>
        )}

        {page === "history" && (
          <div>
            <AnalysisHistory />
          </div>
        )}

        {page === "search" && (
          <div>
            <MoleculeSearch onSelectSmiles={(smi) => { setAnalyzeSmiles(smi); setPage("analyze"); }} />
          </div>
        )}
      </div>
    </div>
  );
}
