import { useMemo, useState } from "react";

type AnalyzeResponse =
  | {
      smiles: string;
      mw: number;
      logp: number;
      hbd: number;
      hba: number;
      tpsa: number;
      num_atoms: number;
    }
  | { error: string };

export default function App() {
  const [smiles, setSmiles] = useState("CCO");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<AnalyzeResponse | null>(null);
  const [err, setErr] = useState<string | null>(null);

  const isError = useMemo(() => result && "error" in result, [result]);

  async function analyze() {
    setLoading(true);
    setErr(null);
    setResult(null);

    try {
      const res = await fetch("/api/analyze", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles }),
      });

      const data = (await res.json()) as AnalyzeResponse;

      if (!res.ok) {
        setErr(`HTTP ${res.status}`);
      }
      setResult(data);
    } catch (e: any) {
      setErr(e?.message ?? "Request failed");
    } finally {
      setLoading(false);
    }
  }

  return (
    <div style={{ maxWidth: 820, margin: "40px auto", padding: 16 }}>
      <h1 style={{ marginBottom: 8 }}>Chembind — SMILES Analyzer</h1>
      <p style={{ marginTop: 0, opacity: 0.75 }}>
        Enter a SMILES string, send it to the FastAPI RDKit backend, and view descriptors.
      </p>

      <div style={{ display: "flex", gap: 12, marginTop: 16 }}>
        <input
          value={smiles}
          onChange={(e) => setSmiles(e.target.value)}
          placeholder="e.g. CCO"
          style={{
            flex: 1,
            padding: "10px 12px",
            borderRadius: 10,
            border: "1px solid rgba(255,255,255,0.15)",
            background: "rgba(255,255,255,0.05)",
            color: "white",
          }}
        />
        <button
          onClick={analyze}
          disabled={loading || !smiles.trim()}
          style={{
            padding: "10px 14px",
            borderRadius: 10,
            border: "1px solid rgba(255,255,255,0.15)",
            background: loading ? "rgba(255,255,255,0.08)" : "rgba(255,255,255,0.12)",
            color: "white",
            cursor: loading ? "not-allowed" : "pointer",
            minWidth: 120,
          }}
        >
          {loading ? "Analyzing…" : "Analyze"}
        </button>
      </div>

      {err && (
        <div
          style={{
            marginTop: 16,
            padding: 12,
            borderRadius: 10,
            background: "rgba(255,0,0,0.12)",
          }}
        >
          <b>Error:</b> {err}
        </div>
      )}

      {result && (
        <div style={{ marginTop: 16 }}>
          <h2 style={{ marginBottom: 8 }}>Result</h2>

          {!isError && "smiles" in result ? (
            <div
              style={{
                display: "grid",
                gridTemplateColumns: "repeat(3, 1fr)",
                gap: 12,
                marginBottom: 16,
              }}
            >
              <Stat label="MW" value={result.mw} />
              <Stat label="logP" value={result.logp} />
              <Stat label="TPSA" value={result.tpsa} />
              <Stat label="HBD" value={result.hbd} />
              <Stat label="HBA" value={result.hba} />
              <Stat label="Atoms" value={result.num_atoms} />
            </div>
          ) : null}

          <pre
            style={{
              padding: 12,
              borderRadius: 12,
              background: "rgba(255,255,255,0.06)",
              overflowX: "auto",
            }}
          >
            {JSON.stringify(result, null, 2)}
          </pre>
        </div>
      )}
    </div>
  );
}

function Stat({ label, value }: { label: string; value: number }) {
  return (
    <div
      style={{
        padding: 12,
        borderRadius: 12,
        background: "rgba(255,255,255,0.06)",
        border: "1px solid rgba(255,255,255,0.10)",
      }}
    >
      <div style={{ fontSize: 12, opacity: 0.75 }}>{label}</div>
      <div style={{ fontSize: 20, fontWeight: 600 }}>{Number(value).toFixed(3)}</div>
    </div>
  );
}

