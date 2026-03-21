import { useState, useEffect } from "react";
import { apiFetch } from "../api/client";

interface Annotation {
  id: string;
  title: string;
  smiles: string;
  notes: { id: string; text: string; atomIndices?: number[] }[];
  created_at: string;
}

export default function AnnotationsTab({
  onSelectSmiles,
}: {
  onSelectSmiles?: (smiles: string) => void;
}) {
  const [annotations, setAnnotations] = useState<Annotation[]>([]);
  const [loading, setLoading] = useState(true);
  const [title, setTitle] = useState("");
  const [smiles, setSmiles] = useState("");
  const [noteText, setNoteText] = useState("");
  const [selectedId, setSelectedId] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);

  const fetchAnnotations = async () => {
    try {
      const resp = await apiFetch<{ annotations: Annotation[] }>("/api/annotations");
      setAnnotations(resp.annotations);
    } catch {
      // Feature may be disabled
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => { fetchAnnotations(); }, []);

  const handleSave = async () => {
    if (!title.trim() || !smiles.trim()) return;
    setError(null);
    try {
      await apiFetch("/api/annotations", {
        method: "POST",
        body: { title: title.trim(), smiles: smiles.trim(), notes: [] },
      });
      setTitle("");
      setSmiles("");
      fetchAnnotations();
    } catch (e: any) {
      setError(e?.message || "Failed to save");
    }
  };

  const handleAddNote = async (annotationId: string) => {
    if (!noteText.trim()) return;
    const annotation = annotations.find((a) => a.id === annotationId);
    if (!annotation) return;

    const newNote = { id: crypto.randomUUID(), text: noteText.trim(), atomIndices: [] };
    const updatedNotes = [...(annotation.notes || []), newNote];

    try {
      await apiFetch(`/api/annotations/${annotationId}`, {
        method: "PUT",
        body: { notes: updatedNotes },
      });
      setNoteText("");
      fetchAnnotations();
    } catch (e: any) {
      setError(e?.message || "Failed to add note");
    }
  };

  const handleShare = async (annotationId: string) => {
    try {
      const resp = await apiFetch<{ shareId: string; url: string }>(
        `/api/annotations/${annotationId}/share`,
        { method: "POST" }
      );
      await navigator.clipboard.writeText(`${window.location.origin}${resp.url}`);
      alert("Share link copied to clipboard!");
    } catch (e: any) {
      setError(e?.message || "Share failed");
    }
  };

  if (loading) return <div style={{ opacity: 0.7 }}>Loading annotations…</div>;

  return (
    <div>
      <div style={{ fontWeight: 600, marginBottom: 8 }}>Annotations</div>

      <div style={{ display: "flex", gap: 8, marginBottom: 12, flexWrap: "wrap" }}>
        <input
          type="text" value={title} onChange={(e) => setTitle(e.target.value)}
          placeholder="Title" style={{ padding: "4px 8px" }}
        />
        <input
          type="text" value={smiles} onChange={(e) => setSmiles(e.target.value)}
          placeholder="SMILES" style={{ padding: "4px 8px", fontFamily: "monospace" }}
        />
        <button onClick={handleSave} disabled={!title.trim() || !smiles.trim()}>
          Save Annotation
        </button>
      </div>

      {error && <div style={{ color: "#c00", fontSize: 13, marginBottom: 8 }}>{error}</div>}

      {annotations.length === 0 && <div style={{ opacity: 0.6, fontSize: 13 }}>No annotations yet.</div>}

      {annotations.map((ann) => (
        <div
          key={ann.id}
          style={{
            border: "1px solid #ddd", borderRadius: 6, padding: 10, marginBottom: 8,
            background: selectedId === ann.id ? "#f0f7ff" : "#fff",
          }}
        >
          <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
            <div>
              <span style={{ fontWeight: 600 }}>{ann.title}</span>
              <span
                style={{ fontFamily: "monospace", fontSize: 12, marginLeft: 8, opacity: 0.7, cursor: "pointer" }}
                onClick={() => onSelectSmiles?.(ann.smiles)}
              >
                {ann.smiles}
              </span>
            </div>
            <div style={{ display: "flex", gap: 4 }}>
              <button onClick={() => setSelectedId(selectedId === ann.id ? null : ann.id)} style={{ fontSize: 11 }}>
                {selectedId === ann.id ? "Collapse" : "Expand"}
              </button>
              <button onClick={() => handleShare(ann.id)} style={{ fontSize: 11 }}>Share</button>
            </div>
          </div>

          {selectedId === ann.id && (
            <div style={{ marginTop: 8, paddingLeft: 8, borderLeft: "2px solid #ddd" }}>
              {(ann.notes || []).map((note, i) => (
                <div key={note.id || i} style={{ fontSize: 13, padding: "2px 0" }}>
                  • {note.text}
                </div>
              ))}

              <div style={{ display: "flex", gap: 4, marginTop: 6 }}>
                <input
                  type="text" value={noteText} onChange={(e) => setNoteText(e.target.value)}
                  placeholder="Add note…" style={{ padding: "3px 6px", fontSize: 12, flex: 1 }}
                  onKeyDown={(e) => e.key === "Enter" && handleAddNote(ann.id)}
                />
                <button onClick={() => handleAddNote(ann.id)} style={{ fontSize: 11 }}>Add</button>
              </div>
            </div>
          )}
        </div>
      ))}
    </div>
  );
}
