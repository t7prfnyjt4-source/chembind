import { useState, useEffect } from "react";

const API_BASE = import.meta.env.VITE_API_BASE || "http://127.0.0.1:8000";

interface SharedAnnotation {
  title: string;
  smiles: string;
  notes: { text: string; atomIndices?: number[] }[];
  viewerState?: Record<string, any>;
  shared_at: string;
}

export default function SharedView({ shareId }: { shareId: string }) {
  const [annotation, setAnnotation] = useState<SharedAnnotation | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    fetch(`${API_BASE}/api/shared/${shareId}`)
      .then((res) => {
        if (!res.ok) throw new Error("Not found");
        return res.json();
      })
      .then((data) => setAnnotation(data))
      .catch((e) => setError(e.message))
      .finally(() => setLoading(false));
  }, [shareId]);

  if (loading) return <div style={{ padding: 24 }}>Loading…</div>;
  if (error) return <div style={{ padding: 24, color: "#c00" }}>Error: {error}</div>;
  if (!annotation) return <div style={{ padding: 24 }}>Not found</div>;

  return (
    <div style={{ padding: 24, fontFamily: "system-ui", maxWidth: 800, margin: "0 auto" }}>
      <h2>{annotation.title}</h2>
      <div style={{ fontFamily: "monospace", fontSize: 14, marginBottom: 12 }}>
        {annotation.smiles}
      </div>

      {annotation.notes && annotation.notes.length > 0 && (
        <div style={{ marginBottom: 12 }}>
          <div style={{ fontWeight: 600, marginBottom: 4 }}>Notes</div>
          {annotation.notes.map((note, i) => (
            <div key={i} style={{ fontSize: 14, padding: "2px 0" }}>• {note.text}</div>
          ))}
        </div>
      )}

      <div style={{ fontSize: 12, opacity: 0.6 }}>
        Shared via ChemBind
      </div>
    </div>
  );
}
