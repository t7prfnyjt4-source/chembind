import { useState } from "react";
import { apiFetch } from "../api/client";

type SearchMode = "similarity" | "substructure";

interface SearchResult {
  doc_id: string;
  smiles: string;
  similarity?: number;
  descriptors: Record<string, number>;
}

interface SearchResponse {
  results: SearchResult[];
  mode: string;
}

export default function MoleculeSearch({
  onSelectSmiles,
}: {
  onSelectSmiles?: (smiles: string) => void;
}) {
  const [mode, setMode] = useState<SearchMode>("similarity");
  const [query, setQuery] = useState("");
  const [results, setResults] = useState<SearchResult[] | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const handleSearch = async () => {
    if (!query.trim()) return;
    setLoading(true);
    setError(null);
    setResults(null);

    try {
      const resp = await apiFetch<SearchResponse>("/api/search", {
        method: "POST",
        body: { query: query.trim(), mode },
      });
      setResults(resp.results);
    } catch (e: any) {
      setError(e?.message || "Search failed");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div>
      <div style={{ fontWeight: 600, marginBottom: 8 }}>Molecule Search</div>

      <div style={{ display: "flex", gap: 8, marginBottom: 12 }}>
        <button
          onClick={() => setMode("similarity")}
          style={{
            fontWeight: mode === "similarity" ? 700 : 400,
            textDecoration: mode === "similarity" ? "underline" : "none",
          }}
        >
          Similar
        </button>
        <button
          onClick={() => setMode("substructure")}
          style={{
            fontWeight: mode === "substructure" ? 700 : 400,
            textDecoration: mode === "substructure" ? "underline" : "none",
          }}
        >
          Substructure
        </button>
      </div>

      <div style={{ display: "flex", gap: 8, marginBottom: 12 }}>
        <input
          type="text"
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder={
            mode === "similarity"
              ? "Enter SMILES (e.g. CCO)"
              : "Enter SMARTS pattern (e.g. [OH])"
          }
          style={{ flex: 1, padding: "6px 10px", fontFamily: "monospace" }}
          onKeyDown={(e) => e.key === "Enter" && handleSearch()}
        />
        <button onClick={handleSearch} disabled={loading || !query.trim()}>
          {loading ? "Searching…" : "Search"}
        </button>
      </div>

      {error && (
        <div style={{ color: "#c00", marginBottom: 8, fontSize: 14 }}>
          {error}
        </div>
      )}

      {results !== null && results.length === 0 && (
        <div style={{ opacity: 0.7, fontSize: 14 }}>No results found.</div>
      )}

      {results && results.length > 0 && (
        <div style={{ display: "flex", flexDirection: "column", gap: 6 }}>
          {results.map((r, i) => (
            <div
              key={r.doc_id || i}
              onClick={() => onSelectSmiles?.(r.smiles)}
              style={{
                padding: "8px 12px",
                border: "1px solid #555",
                borderRadius: 6,
                cursor: onSelectSmiles ? "pointer" : "default",
                fontSize: 14,
              }}
            >
              <div style={{ fontFamily: "monospace", marginBottom: 4 }}>
                {r.smiles}
              </div>
              <div style={{ display: "flex", gap: 12, opacity: 0.8, fontSize: 12 }}>
                {r.similarity !== undefined && (
                  <span>Similarity: {(r.similarity * 100).toFixed(1)}%</span>
                )}
                {r.descriptors?.mw !== undefined && (
                  <span>MW: {r.descriptors.mw.toFixed(1)}</span>
                )}
              </div>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
