import { useState } from "react";
import { GoogleAuthProvider, signInWithPopup, signOut } from "firebase/auth";
import { auth } from "./firebase/config";

export default function App() {
  const [token, setToken] = useState<string>("");
  const [status, setStatus] = useState<string>("");

  const login = async () => {
    setStatus("Logging in...");
    const provider = new GoogleAuthProvider();
    await signInWithPopup(auth, provider);
    const t = await auth.currentUser?.getIdToken(true);
    setToken(t || "");
    setStatus(t ? "Token ready ✅" : "Logged in, but token missing");
  };

  const logout = async () => {
    await signOut(auth);
    setToken("");
    setStatus("Logged out");
  };

  const copy = async () => {
    if (!token) return;
    await navigator.clipboard.writeText(token);
    setStatus("Copied ✅");
  };

  return (
    <div style={{ padding: 40, color: "white", fontFamily: "system-ui" }}>
      <h1>Firebase Token Getter 🔑</h1>

      <div style={{ display: "flex", gap: 12, marginTop: 16 }}>
        <button onClick={login}>Login (Google)</button>
        <button onClick={logout}>Logout</button>
        <button onClick={copy} disabled={!token}>
          Copy token
        </button>
      </div>

      <div style={{ marginTop: 12, opacity: 0.85 }}>{status}</div>

      <textarea
        style={{ width: "100%", height: 240, marginTop: 16 }}
        value={token}
        readOnly
        placeholder="Token will appear here after login..."
      />
    </div>
  );
}
