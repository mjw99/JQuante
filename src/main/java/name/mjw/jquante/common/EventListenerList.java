package name.mjw.jquante.common;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Is a replacement for javax.swing.event.EventListenerList to avoid dependence
 * on swing for core MeTA APIs
 * 
 * @author V. Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class EventListenerList<T> {

	private final HashMap<Class<T>, ArrayList<T>> listenerList;

	/** Creates a new instance of EventListenerList */
	public EventListenerList() {
		listenerList = new HashMap<>();
	}

	/**
	 * Add a new listener
	 * 
	 * @param theClass    the class of listener
	 * @param theListener the listener instance
	 */
	public void add(final Class<T> theClass, final T theListener) {
		if (!listenerList.containsKey(theClass)) {
			listenerList.put(theClass, new ArrayList<T>());
		} // end if

		listenerList.get(theClass).add(theListener);
	}

	/**
	 * Remove a listener
	 * 
	 * @param theClass    the class of listener
	 * @param theListener the listener instance
	 */
	public void remove(final Class<T> theClass, final T theListener) {
		if (!listenerList.containsKey(theClass))
			return;

		listenerList.get(theClass).remove(theListener);
	}

	/**
	 * Get listener list
	 * 
	 * @return the list of listener objects
	 */
	public Object[] getListenerList() {
		final ArrayList<T> listeners = new ArrayList<T>();

		for (final Class<T> listenerClass : listenerList.keySet()) {
			listeners.addAll(listenerList.get(listenerClass));
		}

		return listeners.toArray();
	}

	/**
	 * Get the listener list of a particular class
	 * 
	 * @param theListenerClass the query class
	 * @return the list of listener objects
	 */
	public Object[] getListenerList(final Class<T> theListenerClass) {
		final ArrayList<T> listeners = new ArrayList<T>();

		for (final Class<T> listenerClass : listenerList.keySet()) {
			if (listenerClass.equals(theListenerClass)) {
				listeners.addAll(listenerList.get(listenerClass));
			} // end if
		} // end for

		return listeners.toArray();
	}

	/**
	 * Get the total number of listeners in this list
	 * 
	 * @return the number of listeners in this list
	 */
	public int getListenerCount() {
		int listenerCount = 0;

		for (final Class<T> listenerClass : listenerList.keySet()) {
			listenerCount += listenerList.get(listenerClass).size();
		}

		return listenerCount;
	}

	/**
	 * Get the listener count of a particular class
	 * 
	 * @param theListenerClass the query class
	 * @return number of listeners
	 */
	public int getListenerCount(final Class<T> theListenerClass) {
		int listenerCount = 0;

		for (final Class<T> listenerClass : listenerList.keySet()) {
			if (listenerClass.equals(theListenerClass)) {
				listenerCount += listenerList.get(listenerClass).size();
			}
		}

		return listenerCount;
	}
}
