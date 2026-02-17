import { useState } from "react";
import { GoogleAuthProvider, signInWithPopup, signOut } from "firebase/auth";
import { auth } from "./firebase/config";
import { useAnalyze } from "./hooks/useAnalyze";

function Metric({ label, value }: { label: string; value: string }) {
  return (
    <div style={{ padding: 12, borderRadius: 12, border: "1px solid #e5e7eb" }}>
      <div style={{ fontSize: 12, opacity: 0.6 }}>{label}</div>
      <div style={{ fontSize: 20, fontWeight: 700 }}>{value}</div>
    </div>
  );
}

export default function App() {
  const [smiles, setSmiles] = useState("");
  const { data, loading, error, analyze } = useAnalyze();

  const [authOut, setAuthOut] = useState<any>(null);
  const user = auth.currentUser;

  const onSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    const s = smiles.trim();
    if (s) analyze(s);
  };

  // ---- Auth actions ----
  const login = async () => {
    const provider = new GoogleAuthProvider();
    await signInWithPopup(auth, provider);
    setAuthOut({ loggedIn: true, email: auth.currentUser?.email });
  };

  const logout = async () => {
    await signOut(auth);
    setAuthOut({ loggedOut: true });
  };

  const testMe = async () => {
    const u = auth.currentUser;
    if (!u) return alert("Login first");
    const token = await u.getIdToken();

    const res = await fetch(`${import.meta.env.VITE_API_BASE}/me`, {
      headers: { Authorization: `Bearer ${token}` },
    });

    // Show whatever backend returns
    const text = await res.text();
    try {
      setAuthOut(JSON.parse(text));
    } catch {
      setAuthOut({ raw: text });
    }
  };

  return (
    <div style={{ fontFamily: "system-ui", padding: 24, maxWidth: 900, margin: "0 auto" }}>
      <h1 style={{ fontSize: 34, marginBottom: 8 }}>ChemBind</h1>

      {/* ---- Auth Bar ---- */}
      <div
        style={{
          display: "flex",
          gap: 12,
          alignItems: "center",
          justifyContent: "space-between",
          padding: 12,
          borderRadius: 12,
          border: "1px solid #e5e7eb",
          marginBottom: 16,
        }}
      >
        <div style={{ fontSize: 14, opacity: 0.8 }}>
          {user ? (
            <>
              Signed in as <b>{user.email}</b>
            </>
          ) : (
            <>Not signed in</>
          )}
        </div>

        <div style={{ display: "flex", gap: 10 }}>
          {!user ? (
            <button
              onClick={login}
              style={{ padding: "10px 14px", borderRadius: 10, border: "1px solid #ccc" }}
            >
              Login with Google
            </button>
          ) : (
            <>
              <button
                onClick={testMe}
                style={{ padding: "10px 14px", borderRadius: 10, border: "1px solid #ccc" }}
              >
                Test /me
              </button>
              <button
                onClick={logout}
                style={{ padding: "10px 14px", borderRadius: 10, border: "1px solid #ccc" }}
              >
                Logout
              </button>
            </>
          )}
        </div>
      </div>

      {/* Auth output */}
      {authOut && (
        <pre
          style={{
            marginBottom: 16,
            padding: 12,
            borderRadius: 12,
            border: "1px solid #e5e7eb",
            background: "#0b1220",
            color: "#d1fae5",
            overflowX: "auto",
          }}
        >
          {JSON.stringify(authOut, null, 2)}
        </pre>
      )}

      <div style={{ opacity: 0.7, marginBottom: 16 }}>
        Enter a SMILES string to compute RDKit descriptors.
      </div>

      <form onSubmit={onSubmit} style={{ display: "flex", gap: 12, marginBottom: 16 }}>
        <input
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          placeholder="e.g. CCO, c1ccccc1, CCO.Cl"
          disabled={loading}
          style={{ flex: 1, padding: 12, borderRadius: 10, border: "1px solid #ccc" }}
        />
        <button
          type="submit"
          disabled={loading || !smiles.trim()}
          style={{ padding: "12px 16px", borderRadius: 10, border: "1px solid #ccc" }}
        >
          {loading ? "Analyzing…" : "Analyze"}
        </button>
      </form>

      {error && (
        <div
          style={{
            padding: 12,
            borderRadius: 10,
            border: "1px solid #fca5a5",
            background: "#fef2f2",
            marginBottom: 16,
          }}
        >
          <div style={{ fontWeight: 800 }}>Network / Backend Error</div>
          <div>{error}</div>
        </div>
      )}

      {loading && (
        <div
          style={{
            padding: 12,
            borderRadius: 10,
            border: "1px solid #93c5fd",
            background: "#eff6ff",
            marginBottom: 16,
          }}
        >
          Analyzing molecule…
        </div>
      )}

      {data && !data.valid && (
        <div
          style={{
            padding: 12,
            borderRadius: 10,
            border: "1px solid #fcd34d",
            background: "#fffbeb",
            marginBottom: 16,
          }}
        >
          <div style={{ fontWeight: 800 }}>Invalid SMILES</div>
          <div>{data.error}</div>
          <div style={{ marginTop: 8, fontSize: 12, opacity: 0.6 }}>Request ID: {data.requestId}</div>
        </div>
      )}

      {data && data.valid && data.descriptors && (
        <div style={{ padding: 16, borderRadius: 12, border: "1px solid #e5e7eb" }}>
          <div style={{ marginBottom: 12 }}>
            <div style={{ fontSize: 12, opacity: 0.6 }}>Original SMILES</div>
            <div style={{ fontFamily: "monospace", fontSize: 18 }}>{data.smiles}</div>
          </div>

          {data.canonicalSmiles && data.canonicalSmiles !== data.smiles && (
            <div style={{ marginBottom: 12 }}>
              <div style={{ fontSize: 12, opacity: 0.6 }}>Canonical SMILES</div>
              <div style={{ fontFamily: "monospace", fontSize: 18, color: "#16a34a" }}>
                {data.canonicalSmiles}
              </div>
            </div>
          )}

          <div style={{ display: "grid", gridTemplateColumns: "repeat(2, minmax(0, 1fr))", gap: 12 }}>
            <Metric label="Molecular Weight" value={`${data.descriptors.mw.toFixed(2)} Da`} />
            <Metric label="LogP" value={data.descriptors.logp.toFixed(2)} />
            <Metric label="H-Bond Donors" value={`${data.descriptors.hbd}`} />
            <Metric label="H-Bond Acceptors" value={`${data.descriptors.hba}`} />
            <Metric label="TPSA" value={`${data.descriptors.tpsa.toFixed(2)}`} />
            <Metric label="Atom Count" value={`${data.descriptors.atom_count}`} />
          </div>

          <div style={{ marginTop: 12, fontSize: 12, opacity: 0.6 }}>Request ID: {data.requestId}</div>
        </div>
      )}
    </div>
  );
}

